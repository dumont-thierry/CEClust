#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
scenario_arg <- if (length(args) >= 1L) args[[1L]] else "all"
reps <- if (length(args) >= 2L) as.integer(args[[2L]]) else 5L

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

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
    gauss_univ = list(
      name = "gauss_univ",
      Z = gaussian_mixture_1d(),
      familyType = "gaussUniv",
      r0 = 16L,
      Nshots = 20L,
      Nloop = 160L,
      lambda = 1,
      C = 1
    ),
    gauss_univ_heavy = list(
      name = "gauss_univ_heavy",
      Z = gaussian_mixture_1d(n = 5000),
      familyType = "gaussUniv",
      r0 = 24L,
      Nshots = 50L,
      Nloop = 300L,
      lambda = 1,
      C = 1
    ),
    gauss_univ_single_shot_heavy = list(
      name = "gauss_univ_single_shot_heavy",
      Z = gaussian_mixture_1d(n = 50000),
      familyType = "gaussUniv",
      r0 = 24L,
      Nshots = 1L,
      Nloop = 300L,
      lambda = 1,
      C = 1
    ),
    gauss_vector = list(
      name = "gauss_vector",
      Z = iris[, 1:4],
      familyType = "gaussVector",
      r0 = 8L,
      Nshots = 20L,
      Nloop = 160L,
      lambda = 1,
      C = 10
    ),
    mixed = list(
      name = "mixed",
      Z = simulate_multidim_benchmark_data(n = 500, p_num = 4, p_fac = 2, seed = 11),
      familyType = "gaussAndDiscreteVector",
      r0 = 10L,
      Nshots = 16L,
      Nloop = 120L,
      lambda = 1,
      C = 10
    ),
    stop("Unknown scenario: ", name)
  )
}

run_fit <- function(case, mode, seed) {
  backend_data <- configure_mode(case, mode)

  set.seed(seed)
  gc()
  elapsed <- system.time({
    fit <- CECclassif(
      Z = case$Z,
      lambda = case$lambda,
      C = case$C,
      r0 = case$r0,
      Nshots = case$Nshots,
      Nloop = case$Nloop,
      familyType = case$familyType,
      displayRemainingTime = FALSE,
      backend_data = backend_data
    )
  })[["elapsed"]]

  data.frame(
    scenario = case$name,
    mode = mode,
    seed = seed,
    elapsed = unname(elapsed),
    REO = fit$REO %||% NA_integer_,
    Hphi = fit$Hphi %||% NA_real_,
    stringsAsFactors = FALSE
  )
}

configure_mode <- function(case, mode) {
  if (identical(mode, "current")) {
    CECconfigure_runtime("fast")
    backend_data <- NULL
  } else if (identical(mode, "legacy")) {
    CECconfigure_runtime("base")
    backend_data <- list(optimized = FALSE, raw = case$Z)
  } else if (identical(mode, "fast")) {
    CECconfigure_runtime("fast")
    backend_data <- NULL
  } else {
    stop("Unknown mode: ", mode)
  }
  backend_data
}

warm_up_fit <- function(case, mode, seed) {
  backend_data <- configure_mode(case, mode)
  set.seed(seed)
  gc()
  invisible(
    CECclassif(
      Z = case$Z,
      lambda = case$lambda,
      C = case$C,
      r0 = case$r0,
      Nshots = min(case$Nshots, 2L),
      Nloop = min(case$Nloop, 5L),
      familyType = case$familyType,
      displayRemainingTime = FALSE,
      backend_data = backend_data
    )
  )
}

summarise_timings <- function(raw) {
  do.call(
    rbind,
    lapply(split(raw, list(raw$scenario, raw$mode), drop = TRUE), function(x) {
      data.frame(
        scenario = x$scenario[[1L]],
        mode = x$mode[[1L]],
        reps = nrow(x),
        mean_elapsed = mean(x$elapsed),
        sd_elapsed = if (nrow(x) > 1L) stats::sd(x$elapsed) else NA_real_,
        median_elapsed = median(x$elapsed),
        min_elapsed = min(x$elapsed),
        q25_elapsed = unname(stats::quantile(x$elapsed, 0.25)),
        q75_elapsed = unname(stats::quantile(x$elapsed, 0.75)),
        max_elapsed = max(x$elapsed),
        median_REO = median(x$REO),
        median_Hphi = median(x$Hphi),
        stringsAsFactors = FALSE
      )
    })
  )
}

add_speedups <- function(summary) {
  out <- summary
  out$speedup_vs_legacy <- NA_real_
  for (sc in unique(out$scenario)) {
    legacy_time <- out$mean_elapsed[out$scenario == sc & out$mode == "legacy"]
    if (length(legacy_time) == 1L && is.finite(legacy_time)) {
      idx <- out$scenario == sc
      out$speedup_vs_legacy[idx] <- legacy_time / out$mean_elapsed[idx]
    }
  }
  out
}

modes_for_case <- function(case) {
  if (CECis_fast_family(case$familyType)) {
    c("legacy", "fast")
  } else {
    "current"
  }
}

write_outputs <- function(raw, summary) {
  out_dir <- file.path(pkg_root, "inst", "bench", "out")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  raw_file <- file.path(out_dir, "cecclassif-timings-raw.csv")
  summary_file <- file.path(out_dir, "cecclassif-timings-summary.csv")
  md_file <- file.path(repo_root, "discussions", "accelerer-cecclassif", "benchmark-timings.md")

  write.csv(raw, raw_file, row.names = FALSE)
  write.csv(summary, summary_file, row.names = FALSE)

  lines <- c(
    "# Timings `CECclassif`",
    "",
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("R: ", R.version.string),
    paste0("Repetitions: ", reps),
    "",
    "Commande:",
    "",
    "```sh",
    paste0("Rscript CEClust-main/inst/bench/benchmark-cecclassif.R ", scenario_arg, " ", reps),
    "```",
    "",
    "Metrique principale: `mean_elapsed`.",
    "",
    "Note: ce script utilise `pkgload::load_all()`. Pour juger les optimisations",
    "C++ compilees, utiliser `benchmark-cecclassif-installed.R`, qui installe",
    "une copie temporaire du package avant mesure.",
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

scenarios <- if (identical(scenario_arg, "all")) {
  c("gauss_univ", "gauss_vector", "mixed")
} else {
  strsplit(scenario_arg, ",", fixed = TRUE)[[1L]]
}

results <- list()
for (sc in scenarios) {
  case <- make_case(sc)
  modes <- modes_for_case(case)
  for (mode in modes) {
    warm_up_fit(case, mode, seed = 20260525L)
  }
  for (rep_idx in seq_len(reps)) {
    seed <- 20260526L + rep_idx
    for (mode in modes) {
      row <- run_fit(case, mode, seed)
      print(row)
      results[[length(results) + 1L]] <- row
    }
  }
}

raw <- do.call(rbind, results)
summary <- add_speedups(summarise_timings(raw))
md_file <- write_outputs(raw, summary)

cat("\nTiming summary written to ", md_file, "\n", sep = "")
