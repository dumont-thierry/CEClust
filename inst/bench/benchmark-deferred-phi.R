#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
scenario_arg <- if (length(args) >= 1L) args[[1L]] else "mixed"
reps <- if (length(args) >= 2L) as.integer(args[[2L]]) else 20L

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

make_case <- function(name) {
  switch(
    name,
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
    stop("Unknown scenario: ", name, ". This benchmark is currently meaningful only for 'mixed'.")
  )
}

set_binding <- function(env, name, value) {
  was_locked <- bindingIsLocked(name, env)
  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }
}

install_dense_phi_wrapper <- function() {
  ns <- environment(CECclassif)
  original <- get("CECoptPhi_backend", envir = ns, inherits = FALSE)
  wrapper <- function(Z, param, lambda = 1, backend_data = NULL) {
    out <- original(Z = Z, param = param, lambda = lambda, backend_data = backend_data)
    if (is.null(out$phi) && !is.null(out$clusters)) {
      out$phi <- CECclusters_to_phi(out$clusters, length(param$nu))
    }
    out
  }
  set_binding(ns, "CECoptPhi_backend", wrapper)
  function() set_binding(ns, "CECoptPhi_backend", original)
}

run_fit <- function(case, mode, seed) {
  CECconfigure_runtime("fast")
  backend_data <- CECprepare_backend_data(case$Z, familyType = case$familyType)
  restore <- NULL
  if (identical(mode, "dense_phi")) {
    restore <- install_dense_phi_wrapper()
    on.exit(restore(), add = TRUE)
  }

  set.seed(seed)
  gc()
  fit <- NULL
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

summarise_timings <- function(raw) {
  do.call(
    rbind,
    lapply(split(raw, raw$mode, drop = TRUE), function(x) {
      data.frame(
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

add_speedup <- function(summary) {
  dense_time <- summary$mean_elapsed[summary$mode == "dense_phi"]
  if (length(dense_time) == 1L && is.finite(dense_time)) {
    summary$speedup_vs_dense_phi <- dense_time / summary$mean_elapsed
  } else {
    summary$speedup_vs_dense_phi <- NA_real_
  }
  summary
}

write_outputs <- function(raw, summary) {
  out_dir <- file.path(pkg_root, "inst", "bench", "out")
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  scenario_slug <- gsub("[^A-Za-z0-9_-]+", "-", scenario_arg)
  raw_file <- file.path(out_dir, paste0("cecclassif-deferred-phi-", scenario_slug, "-raw.csv"))
  summary_file <- file.path(out_dir, paste0("cecclassif-deferred-phi-", scenario_slug, "-summary.csv"))
  md_file <- file.path(repo_root, "discussions", "accelerer-cecclassif", paste0("deferred-phi-", scenario_slug, ".md"))
  latest_md_file <- file.path(repo_root, "discussions", "accelerer-cecclassif", "deferred-phi.md")

  write.csv(raw, raw_file, row.names = FALSE)
  write.csv(summary, summary_file, row.names = FALSE)

  lines <- c(
    "# Benchmark `phi` differe",
    "",
    paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
    paste0("R: ", R.version.string),
    "",
    "Commande:",
    "",
    "```sh",
    paste0("Rscript CEClust-main/inst/bench/benchmark-deferred-phi.R ", scenario_arg, " ", reps),
    "```",
    "",
    "Metrique principale: `mean_elapsed`.",
    "",
    "Ce benchmark est borne au scenario `mixed`, car le report de `phi` dense est maintenant limite au chemin rapide mixte.",
    "",
    "Modes:",
    "",
    "- `deferred_phi`: chemin courant, `phi` public reconstruit seulement a la fin;",
    "- `dense_phi`: wrapper de controle qui reconstruit `phi` a chaque iteration optimisee.",
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
  writeLines(lines, latest_md_file)
  md_file
}

case <- make_case(scenario_arg)
rows <- list()
for (rep_idx in seq_len(reps)) {
  seed <- 20260526L + rep_idx
  for (mode in c("deferred_phi", "dense_phi")) {
    row <- run_fit(case, mode, seed)
    print(row)
    rows[[length(rows) + 1L]] <- row
  }
}

raw <- do.call(rbind, rows)
summary <- add_speedup(summarise_timings(raw))
md_file <- write_outputs(raw, summary)
cat("\nDeferred phi summary written to ", md_file, "\n", sep = "")
