#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
scenario_arg <- if (length(args) >= 1L) args[[1L]] else "all"
reps <- if (length(args) >= 2L) as.integer(args[[2L]]) else 5L

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
tmp_root <- tempfile("cec-installed-bench-")
src_root <- file.path(tmp_root, "src")
install_lib <- file.path(tmp_root, "lib")
dir.create(src_root, recursive = TRUE)
dir.create(install_lib, recursive = TRUE)

copy_ok <- file.copy(pkg_root, tmp_root, recursive = TRUE)
if (!copy_ok) {
  stop("Could not copy package to temporary benchmark directory.")
}
pkg_copy <- file.path(tmp_root, basename(pkg_root))
src_dir <- file.path(pkg_copy, "src")
if (dir.exists(src_dir)) {
  unlink(
    list.files(src_dir, pattern = "\\.(o|so|dll)$", full.names = TRUE),
    force = TRUE
  )
}

lib_paths <- c(install_lib)
if (dir.exists(local_lib)) {
  lib_paths <- c(lib_paths, normalizePath(local_lib, mustWork = TRUE))
}
lib_paths <- c(lib_paths, .libPaths())
env <- c(
  paste0("R_LIBS=", paste(lib_paths, collapse = .Platform$path.sep)),
  paste0("R_LIBS_USER=", if (dir.exists(local_lib)) normalizePath(local_lib, mustWork = TRUE) else "")
)

install_log <- system2(
  file.path(R.home("bin"), "R"),
  c("CMD", "INSTALL", "--preclean", "--no-docs", "--no-help", "--no-html", "-l", install_lib, pkg_copy),
  stdout = TRUE,
  stderr = TRUE,
  env = env
)
install_status <- attr(install_log, "status")
if (!is.null(install_status) && install_status != 0L) {
  writeLines(install_log)
  stop("R CMD INSTALL failed.")
}

runner <- file.path(tmp_root, "run-installed-benchmark.R")
writeLines(c(
  "args <- commandArgs(trailingOnly = TRUE)",
  "scenario_arg <- args[[1L]]",
  "reps <- as.integer(args[[2L]])",
  "install_lib <- args[[3L]]",
  "local_lib <- args[[4L]]",
  ".libPaths(c(install_lib, local_lib, .libPaths()))",
  "library(CEClust)",
  "`%||%` <- function(x, y) if (is.null(x)) y else x",
  "gaussian_mixture_1d <- function(n = 800L, seed = 13L) {",
  "  set.seed(seed)",
  "  mu <- 2 * sqrt(3) * c(-3, -1.75, -0.75, 0, 0.75, 1.75, 3)",
  "  sd_vec <- c(0.2, 1, 0.5, 1, 0.5, 1, 0.2)",
  "  pi_vec <- c(0.5, 2, 0.5, 1, 0.5, 2, 0.5)",
  "  pi_vec <- pi_vec / sum(pi_vec)",
  "  component <- sample.int(length(mu), size = n, replace = TRUE, prob = pi_vec)",
  "  rnorm(n, mean = mu[component], sd = sd_vec[component])",
  "}",
  "simulate_multidim <- function(...) {",
  "  get('simulate_multidim_benchmark_data', envir = asNamespace('CEClust'))(...)",
  "}",
  "make_case <- function(name) {",
  "  switch(",
  "    name,",
  "    gauss_univ = list(name = 'gauss_univ', Z = gaussian_mixture_1d(), familyType = 'gaussUniv', r0 = 16L, Nshots = 20L, Nloop = 160L, lambda = 1, C = 1),",
  "    gauss_univ_heavy = list(name = 'gauss_univ_heavy', Z = gaussian_mixture_1d(n = 5000L), familyType = 'gaussUniv', r0 = 24L, Nshots = 50L, Nloop = 300L, lambda = 1, C = 1),",
  "    gauss_univ_single_shot_heavy = list(name = 'gauss_univ_single_shot_heavy', Z = gaussian_mixture_1d(n = 50000L), familyType = 'gaussUniv', r0 = 24L, Nshots = 1L, Nloop = 300L, lambda = 1, C = 1),",
  "    gauss_vector = list(name = 'gauss_vector', Z = iris[, 1:4], familyType = 'gaussVector', r0 = 8L, Nshots = 20L, Nloop = 160L, lambda = 1, C = 10),",
  "    mixed = list(name = 'mixed', Z = simulate_multidim(n = 500, p_num = 4, p_fac = 2, seed = 11), familyType = 'gaussAndDiscreteVector', r0 = 10L, Nshots = 16L, Nloop = 120L, lambda = 1, C = 10),",
  "    stop('Unknown scenario: ', name)",
  "  )",
  "}",
  "configure_mode <- function(case, mode) {",
  "  if (identical(mode, 'current')) {",
  "    CECconfigure_runtime('fast')",
  "    NULL",
  "  } else if (identical(mode, 'legacy')) {",
  "    CECconfigure_runtime('base')",
  "    list(optimized = FALSE, raw = case$Z)",
  "  } else if (identical(mode, 'fast')) {",
  "    CECconfigure_runtime('fast')",
  "    NULL",
  "  } else {",
  "    stop('Unknown mode: ', mode)",
  "  }",
  "}",
  "modes_for_case <- function(case) {",
  "  if (case$familyType %in% c('gaussVector', 'discreteVector', 'gaussAndDiscreteVector')) c('legacy', 'fast') else 'current'",
  "}",
  "run_fit <- function(case, mode, seed) {",
  "  backend_data <- configure_mode(case, mode)",
  "  set.seed(seed)",
  "  gc()",
  "  elapsed <- system.time({",
  "    fit <- CECclassif(Z = case$Z, lambda = case$lambda, C = case$C, r0 = case$r0, Nshots = case$Nshots, Nloop = case$Nloop, familyType = case$familyType, displayRemainingTime = FALSE, backend_data = backend_data)",
  "  })[['elapsed']]",
  "  data.frame(scenario = case$name, mode = mode, seed = seed, elapsed = unname(elapsed), REO = fit$REO %||% NA_integer_, Hphi = fit$Hphi %||% NA_real_, stringsAsFactors = FALSE)",
  "}",
  "summarise_timings <- function(raw) {",
  "  do.call(rbind, lapply(split(raw, list(raw$scenario, raw$mode), drop = TRUE), function(x) {",
  "    data.frame(scenario = x$scenario[[1L]], mode = x$mode[[1L]], reps = nrow(x), mean_elapsed = mean(x$elapsed), sd_elapsed = if (nrow(x) > 1L) stats::sd(x$elapsed) else NA_real_, median_elapsed = stats::median(x$elapsed), min_elapsed = min(x$elapsed), q25_elapsed = unname(stats::quantile(x$elapsed, 0.25)), q75_elapsed = unname(stats::quantile(x$elapsed, 0.75)), max_elapsed = max(x$elapsed), median_REO = stats::median(x$REO), median_Hphi = stats::median(x$Hphi), stringsAsFactors = FALSE)",
  "  }))",
  "}",
  "add_speedups <- function(summary) {",
  "  out <- summary",
  "  out$speedup_vs_legacy <- NA_real_",
  "  for (sc in unique(out$scenario)) {",
  "    legacy_time <- out$mean_elapsed[out$scenario == sc & out$mode == 'legacy']",
  "    if (length(legacy_time) == 1L && is.finite(legacy_time)) {",
  "      idx <- out$scenario == sc",
  "      out$speedup_vs_legacy[idx] <- legacy_time / out$mean_elapsed[idx]",
  "    }",
  "  }",
  "  out",
  "}",
  "scenarios <- if (identical(scenario_arg, 'all')) c('gauss_univ', 'gauss_vector', 'mixed') else strsplit(scenario_arg, ',', fixed = TRUE)[[1L]]",
  "rows <- list()",
  "for (scenario in scenarios) {",
  "  case <- make_case(scenario)",
  "  for (mode in modes_for_case(case)) {",
  "    invisible(run_fit(case, mode, 20260526L))",
  "    for (rep_idx in seq_len(reps)) {",
  "      row <- run_fit(case, mode, 20260526L + rep_idx)",
  "      print(row)",
  "      rows[[length(rows) + 1L]] <- row",
  "    }",
  "  }",
  "}",
  "raw <- do.call(rbind, rows)",
  "summary <- add_speedups(summarise_timings(raw))",
  "saveRDS(list(raw = raw, summary = summary), file = file.path(tempdir(), 'cec-installed-benchmark.rds'))",
  "print(summary, row.names = FALSE)"
), runner)

runner_out <- system2(
  file.path(R.home("bin"), "Rscript"),
  c(runner, scenario_arg, as.character(reps), install_lib, if (dir.exists(local_lib)) normalizePath(local_lib, mustWork = TRUE) else ""),
  stdout = TRUE,
  stderr = TRUE,
  env = env
)
runner_status <- attr(runner_out, "status")
if (!is.null(runner_status) && runner_status != 0L) {
  writeLines(runner_out)
  stop("Installed benchmark failed.")
}

summary_start <- max(grep("^\\s*scenario\\s+mode\\s+reps", runner_out), na.rm = TRUE)
summary_lines <- if (is.finite(summary_start)) runner_out[summary_start:length(runner_out)] else runner_out

out_dir <- file.path(pkg_root, "inst", "bench", "out")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
md_file <- file.path(repo_root, "discussions", "accelerer-cecclassif", "benchmark-installed-timings.md")
install_log_file <- file.path(out_dir, "cecclassif-installed-install.log")
writeLines(install_log, install_log_file)

lines <- c(
  "# Timings installes `CECclassif`",
  "",
  paste0("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("R: ", R.version.string),
  paste0("Repetitions: ", reps),
  "",
  "Commande:",
  "",
  "```sh",
  paste0("Rscript CEClust-main/inst/bench/benchmark-cecclassif-installed.R ", scenario_arg, " ", reps),
  "```",
  "",
  "Metrique principale: `mean_elapsed`.",
  "",
  "Ce benchmark installe une copie temporaire du package avec `R CMD INSTALL`",
  "avant de mesurer. La copie temporaire est nettoyee avant installation pour",
  "eviter de reutiliser des objets `pkgload` compiles en `-O0`. Il est le",
  "repere adapte aux optimisations C++.",
  "",
  "## Resume",
  "",
  "```",
  summary_lines,
  "```",
  "",
  "## Compilation",
  "",
  "```text",
  grep("g\\+\\+|clang\\+\\+|cec_fast.cpp|using C\\+\\+ compiler", install_log, value = TRUE),
  "```",
  "",
  paste0("- install log: `", normalizePath(install_log_file, mustWork = FALSE), "`")
)
writeLines(lines, md_file)
writeLines(runner_out)
cat("\nInstalled timing summary written to ", md_file, "\n", sep = "")
