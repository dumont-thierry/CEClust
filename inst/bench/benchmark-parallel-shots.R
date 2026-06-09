#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
scenario_arg <- if (length(args) >= 1L) args[[1L]] else "gauss_univ_heavy"
reps <- if (length(args) >= 2L) as.integer(args[[2L]]) else 10L
n_cores <- if (length(args) >= 3L) as.integer(args[[3L]]) else 2L

find_pkg_root <- function() {
  candidates <- c(
    normalizePath(getwd(), mustWork = FALSE),
    normalizePath(file.path(getwd(), "CEClust-main"), mustWork = FALSE),
    normalizePath(file.path(getwd(), "..", ".."), mustWork = FALSE)
  )
  for (candidate in unique(candidates)) {
    if (file.exists(file.path(candidate, "DESCRIPTION")) &&
        basename(candidate) == "CEClust-main") {
      return(candidate)
    }
  }
  stop("Could not locate CEClust-main. Run this script from the repository root or package root.")
}

pkg_root <- find_pkg_root()
repo_root <- dirname(pkg_root)
local_lib <- file.path(repo_root, ".Rlib")
if (dir.exists(local_lib)) {
  .libPaths(c(normalizePath(local_lib, mustWork = TRUE), .libPaths()))
  Sys.setenv(R_LIBS = paste(.libPaths(), collapse = .Platform$path.sep))
  Sys.setenv(R_LIBS_USER = normalizePath(local_lib, mustWork = TRUE))
}

if (!requireNamespace("pkgload", quietly = TRUE)) {
  stop("Package 'pkgload' is required. Install it in .Rlib or the active R library.")
}
pkgload::load_all(pkg_root, quiet = TRUE)

gaussian_mixture_1d <- function(n = 800L, seed = 13L) {
  set.seed(seed)
  mu <- 2 * sqrt(3) * c(-3, -1.75, -0.75, 0, 0.75, 1.75, 3)
  sd_vec <- c(0.2, 1, 0.5, 1, 0.5, 1, 0.2)
  pi_vec <- c(0.5, 2, 0.5, 1, 0.5, 2, 0.5)
  pi_vec <- pi_vec / sum(pi_vec)
  component <- sample.int(length(mu), size = n, replace = TRUE, prob = pi_vec)
  rnorm(n, mean = mu[component], sd = sd_vec[component])
}

make_case <- function(name) {
  switch(
    name,
    gauss_univ_heavy = list(
      name = "gauss_univ_heavy",
      Z = gaussian_mixture_1d(n = 5000),
      familyType = "gaussUniv",
      r0 = 24L,
      Nshots = 50L,
      Nloop = 300L,
      lambda = 1,
      C = 1,
      Cquali = Inf
    ),
    gauss_univ = list(
      name = "gauss_univ",
      Z = gaussian_mixture_1d(n = 800),
      familyType = "gaussUniv",
      r0 = 16L,
      Nshots = 20L,
      Nloop = 160L,
      lambda = 1,
      C = 1,
      Cquali = Inf
    ),
    stop("Unknown scenario: ", name)
  )
}

time_expr <- function(expr) {
  elapsed <- system.time(force(expr))[["elapsed"]]
  unname(elapsed)
}

run_current <- function(case, seed) {
  CECconfigure_runtime("fast")
  set.seed(seed)
  backend_data <- CECprepare_backend_data(case$Z, familyType = case$familyType)
  fit <- NULL
  elapsed <- time_expr({
    fit <- CECclassif(
      Z = case$Z,
      lambda = case$lambda,
      C = case$C,
      Cquali = case$Cquali,
      r0 = case$r0,
      Nshots = case$Nshots,
      Nloop = case$Nloop,
      familyType = case$familyType,
      backend_data = backend_data
    )
  })
  list(elapsed = elapsed, fit = fit)
}

run_independent_sequential <- function(case, seed) {
  CECconfigure_runtime("fast")
  set.seed(seed)
  backend_data <- CECprepare_backend_data(case$Z, familyType = case$familyType)
  fit <- NULL
  elapsed <- time_expr({
    fit <- CECclassif(
      Z = case$Z,
      lambda = case$lambda,
      C = case$C,
      Cquali = case$Cquali,
      r0 = case$r0,
      Nshots = case$Nshots,
      Nloop = case$Nloop,
      familyType = case$familyType,
      backend_data = backend_data,
      shot_strategy = "independent",
      n_cores = 1L
    )
  })
  list(elapsed = elapsed, fit = fit)
}

run_independent_parallel <- function(case, seed, n_cores) {
  CECconfigure_runtime("fast")
  set.seed(seed)
  backend_data <- CECprepare_backend_data(case$Z, familyType = case$familyType)
  n_cores <- max(1L, as.integer(n_cores))
  fit <- NULL
  elapsed <- time_expr({
    fit <- CECclassif(
      Z = case$Z,
      lambda = case$lambda,
      C = case$C,
      Cquali = case$Cquali,
      r0 = case$r0,
      Nshots = case$Nshots,
      Nloop = case$Nloop,
      familyType = case$familyType,
      backend_data = backend_data,
      shot_strategy = "independent",
      n_cores = n_cores
    )
  })
  list(elapsed = elapsed, fit = fit)
}

result_row <- function(case, mode, seed, elapsed, fit, n_cores_used) {
  data.frame(
    scenario = case$name,
    mode = mode,
    seed = seed,
    n_cores = n_cores_used,
    elapsed = elapsed,
    REO = if (!is.null(fit)) fit$REO else NA_integer_,
    Hphi = if (!is.null(fit)) fit$Hphi else NA_real_,
    stringsAsFactors = FALSE
  )
}

summarise_timings <- function(raw) {
  do.call(
    rbind,
    lapply(split(raw, raw$mode, drop = TRUE), function(x) {
      data.frame(
        mode = x$mode[[1L]],
        reps = nrow(x),
        n_cores = x$n_cores[[1L]],
        mean_elapsed = mean(x$elapsed),
        sd_elapsed = if (nrow(x) > 1L) stats::sd(x$elapsed) else NA_real_,
        median_elapsed = median(x$elapsed),
        min_elapsed = min(x$elapsed),
        q25_elapsed = unname(stats::quantile(x$elapsed, 0.25)),
        q75_elapsed = unname(stats::quantile(x$elapsed, 0.75)),
        max_elapsed = max(x$elapsed),
        median_REO = median(x$REO, na.rm = TRUE),
        median_Hphi = median(x$Hphi, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
}

write_outputs <- function(raw, summary) {
  out_dir <- file.path(pkg_root, "inst", "bench", "out")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  raw_file <- file.path(out_dir, "cecclassif-parallel-shots-raw.csv")
  summary_file <- file.path(out_dir, "cecclassif-parallel-shots-summary.csv")
  md_file <- file.path(repo_root, "discussions", "accelerer-cecclassif", "parallel-shots.md")

  write.csv(raw, raw_file, row.names = FALSE)
  write.csv(summary, summary_file, row.names = FALSE)

  lines <- c(
    "# Benchmark parallelisation des shots",
    "",
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("R: ", R.version.string),
    "",
    "Commande:",
    "",
    "```sh",
    paste0("Rscript CEClust-main/inst/bench/benchmark-parallel-shots.R ", scenario_arg, " ", reps, " ", n_cores),
    "```",
    "",
    "Metrique principale: `mean_elapsed`.",
    "",
    "Important: les modes `independent_*` utilisent `CECclassif(..., shot_strategy = \"independent\")`.",
    "Ils lancent des shots `phi0 = NULL` independants et ne reproduisent donc pas exactement la strategie actuelle, qui reutilise parfois `PHI` ou `bestClassif`.",
    "",
    "## Usage recommande",
    "",
    "- garder `shot_strategy = \"current\"` pour reproduire le comportement historique;",
    "- utiliser `shot_strategy = \"independent\"` quand `Nshots` est assez grand et que les shots independants sont acceptables;",
    "- augmenter `n_cores` seulement pour ce mode independant; le mode historique reste sequentiel;",
    "- comparer les gains avec `mean_elapsed`, pas avec une seule execution isolee.",
    "",
    "## Resume",
    "",
    "```",
    capture.output(print(summary, row.names = FALSE)),
    "```",
    "",
    paste0("- raw csv: `", normalizePath(raw_file, mustWork = FALSE), "`"),
    paste0("- summary csv: `", normalizePath(summary_file, mustWork = FALSE), "`")
  )
  writeLines(lines, md_file)
  md_file
}

case <- make_case(scenario_arg)
rows <- list()
for (rep_idx in seq_len(reps)) {
  seed <- 20260526L + rep_idx

  current <- run_current(case, seed)
  row <- result_row(case, "current", seed, current$elapsed, current$fit, 1L)
  print(row)
  rows[[length(rows) + 1L]] <- row

  independent_seq <- run_independent_sequential(case, seed)
  row <- result_row(case, "independent_sequential", seed, independent_seq$elapsed, independent_seq$fit, 1L)
  print(row)
  rows[[length(rows) + 1L]] <- row

  independent_par <- run_independent_parallel(case, seed, n_cores)
  row <- result_row(case, "independent_parallel", seed, independent_par$elapsed, independent_par$fit, n_cores)
  print(row)
  rows[[length(rows) + 1L]] <- row
}

raw <- do.call(rbind, rows)
summary <- summarise_timings(raw)
md_file <- write_outputs(raw, summary)
cat("\nParallel shots summary written to ", md_file, "\n", sep = "")
