bands_progress <- function(store, phase, details = "") {
  value <- list(phase = phase, details = details, updated_at = format(Sys.time(), usetz = TRUE))
  bands_save(value, file.path(store$path, "progress.rds"))
  html <- c("<!doctype html><html lang='en'><meta charset='utf-8'>",
    "<meta http-equiv='refresh' content='20'><title>CEClust progress</title>",
    "<style>body{font:18px system-ui;max-width:900px;margin:50px auto;line-height:1.6}</style>",
    "<h1>CEClust</h1>", paste0("<p>", phase, "</p><p>", details, "</p>"),
    paste0("<p>", value$updated_at, "</p></html>"))
  writeLines(html, file.path(store$path, "SUIVI.html"), useBytes = TRUE)
  invisible(value)
}

cec_bands <- function(X, L = 20L, B = 20L, workers = 1L, seed = 1L,
    bootstrap_indices = NULL, pilot_indices = NULL, control = bands_control(),
    penalty = c("none", "bic"), directory = tempfile("bands-"), resume = TRUE) {
  started <- proc.time()[3L]; rng <- bands_rng_save(); on.exit(bands_rng_restore(rng))
  X <- bands_data(X); penalty <- match.arg(penalty)
  L <- bands_integer(L, "L", 2L); B <- bands_integer(B, "B"); workers <- bands_integer(workers, "workers")
  seed <- bands_integer(seed, "seed", 0L)
  if (!is.logical(resume) || length(resume) != 1L || is.na(resume)) stop("resume must be TRUE or FALSE.")
  # Re-validate a possibly edited control list through the public constructor.
  control <- do.call(bands_control, control)
  store <- bands_store(directory); runtime <- bands_runtime(); n <- nrow(X)
  if (is.null(pilot_indices)) {
    set.seed(bands_seed(seed, "pilot_indices")); pilot_indices <- sample.int(n, n, replace = TRUE)
  }
  bands_check_indices(pilot_indices, n); pilot_indices <- as.integer(pilot_indices)
  jobs <- bands_jobs(n, L, B, seed, bootstrap_indices)
  signature <- digest::digest(list(X, pilot_indices, jobs, control, penalty, runtime))
  manifest_path <- file.path(store$path, "manifest.rds")
  if (file.exists(manifest_path)) {
    previous <- readRDS(manifest_path)
    if (!resume) stop("Archive already exists; choose a new directory or resume=TRUE.")
    if (!identical(previous$signature, signature)) stop("Incompatible data, settings, seeds or runtime for resume.")
  } else {
    contents <- setdiff(list.files(store$path, all.files = TRUE), c(".", "..", "labels"))
    if (length(contents)) stop("New run requires an empty archive directory.")
    bands_save(list(signature = signature, runtime = runtime, control = control,
      penalty = penalty, seed = seed, jobs = jobs, pilot_indices = pilot_indices,
      data_key = digest::digest(X), created_at = format(Sys.time(), usetz = TRUE)), manifest_path)
    bands_save(X, file.path(store$path, "data.rds"))
  }
  final_path <- file.path(store$path, "result.rds")
  if (file.exists(final_path)) {
    result <- readRDS(final_path)
    if (!identical(result$signature, signature)) stop("Result signature mismatch.")
    result$store <- store
    bands_publish_result(result)
    return(result)
  }
  bands_progress(store, "pilot", "Independent bootstrap pilot and paired CEM descents.")
  pilot_path <- file.path(store$path, "pilot.rds")
  if (file.exists(pilot_path)) pilot <- readRDS(pilot_path)
  else {
    pilot <- bands_pilot(X, pilot_indices, bands_seed(seed, "pilot_search"), control, penalty, store)
    bands_save(pilot, pilot_path)
  }
  anchors <- bands_choose_anchors(pilot, control)
  bands_save(anchors, file.path(store$path, "anchors.rds"))
  bands_progress(store, "banks", paste(L, "original banks and", B,
    "bootstrap banks, each with its own internal bootstrap."))
  bank_started <- proc.time()[3L]
  runs <- bands_run_jobs(jobs, X, pilot, anchors, control, penalty, store, runtime, workers)
  bank_seconds <- unname(proc.time()[3L] - bank_started)
  bands_progress(store, "reference", "Pooled original candidates and bootstrap projections; adaptive common diagram.")
  reference_path <- file.path(store$path, "reference.rds")
  if (file.exists(reference_path)) reference <- readRDS(reference_path)
  else {
    reference <- bands_reference(X, runs, anchors, control, penalty, store)
    bands_save(reference, reference_path)
  }
  diagnostic_signature <- digest::digest(list(signature, reference$diagram$tiles,
    names(reference$candidates), vapply(runs, `[[`, "", "signature")))
  bands_progress(store, "diagnostics", "Independent original-original and original-bootstrap agreements on the common diagram.")
  diagnostics <- bands_diagnose(X, runs, reference, L, B, control, penalty, store, diagnostic_signature)
  limits <- unique(c(pilot$diagram$limits, reference$limits,
    if (isTRUE(attr(anchors, "limited"))) "pilot_anchor_budget",
    if (length(pilot$skipped_orders)) "pilot_initial_orders_exceed_sample_capacity",
    unlist(lapply(runs, `[[`, "limits"), use.names = FALSE)))
  result <- structure(list(X = X, L = L, B = B, workers = workers, seed = seed,
    control = control, penalty = penalty, store = store, runtime = runtime, signature = signature,
    pilot = pilot, anchors = anchors, runs = runs, reference = reference,
    tiles = diagnostics$tiles, points = diagnostics$points, limits = limits,
    cache_statistics = diagnostics$cache,
    timings = list(pilot = pilot$seconds, banks_this_call = bank_seconds,
      bank_execution_sum = sum(vapply(runs, `[[`, 0, "seconds")),
      reference = reference$seconds, diagnostics_this_call = diagnostics$seconds,
      total_this_call = unname(proc.time()[3L] - started)),
    scope = "Shared pilot; symmetric estimator on each target sample. Regional witnesses unless diag_mode='all'. Finite-bank exploration."),
    class = "cec_bands")
  bands_save(result, final_path)
  bands_publish_result(result)
  result
}

bands_publish_result <- function(result) {
  store <- result$store
  utils::write.csv(result$tiles, file.path(store$path, "tiles.csv"), row.names = FALSE)
  utils::write.csv(result$points, file.path(store$path, "points.csv"), row.names = FALSE)
  utils::write.csv(bands_partitions(result), file.path(store$path, "partitions.csv"), row.names = FALSE)
  bands_progress(store, "complete", paste(nrow(result$points), "points measured;",
    length(unique(result$tiles$key[result$tiles$keep])), "distinct retained partitions."))
  invisible(result)
}

bands_at <- function(fit, C, lambda) {
  if (!inherits(fit, "cec_bands")) stop("Supply a cec_bands result.")
  if (length(C) != 1L || !is.finite(C) || C <= 0 || length(lambda) != 1L ||
      !is.finite(lambda) || lambda <= 0) stop("C and lambda must be positive finite scalars.")
  cache <- bands_cache(fit$control$cache_mb * 1024^2)
  bands_evaluate_point(fit$X, fit$runs, fit$reference$candidates, C, lambda, fit$L, fit$B,
    fit$control, fit$penalty, fit$store, cache)
}

bands_partitions <- function(fit) {
  tiles <- fit$tiles; keys <- unique(tiles$key[tiles$keep]); rows <- list()
  for (key in keys) {
    ids <- which(tiles$key == key & tiles$keep)
    witness <- tiles$witness[ids[1L]]; point <- fit$points[match(witness, fit$points$id), ]
    candidate <- fit$reference$candidates[[key]]
    labels <- bands_load_labels(candidate, fit$store)
    rows[[length(rows) + 1L]] <- data.frame(key = key, K = candidate$K,
      sizes = paste(sort(tabulate(labels), decreasing = TRUE), collapse = "/"),
      witness = witness, C = point$C, lambda = point$lambda,
      stab_algo = point$stab_algo, stab_boot = point$stab_boot, RSI = point$RSI,
      witness_saturation = tiles$saturation[match(witness, tiles$id)],
      min_retained_saturation = min(tiles$saturation[ids]),
      max_retained_saturation = max(tiles$saturation[ids]), cells = length(ids))
  }
  if (!length(rows)) return(data.frame(name = character(), key = character(), K = integer(),
    sizes = character(), witness = integer(), C = numeric(), lambda = numeric(),
    stab_algo = numeric(), stab_boot = numeric(), RSI = numeric(), witness_saturation = numeric(),
    min_retained_saturation = numeric(), max_retained_saturation = numeric(), cells = integer()))
  result <- do.call(rbind, rows); result <- result[order(result$K, result$key), ]
  result$name <- paste0("P", seq_len(nrow(result))); rownames(result) <- NULL
  result[, c("name", setdiff(names(result), "name"))]
}

bands_plot <- function(fit, filtered = TRUE, main = NULL) {
  if (!inherits(fit, "cec_bands")) stop("Supply a cec_bands result.")
  tiles <- fit$tiles; K <- sort(unique(tiles$K))
  colors <- grDevices::hcl.colors(length(K), "Dark 3")
  fill <- colors[match(tiles$K, K)]
  if (filtered) fill[!tiles$keep] <- "#E1E4E8"
  if (is.null(main)) main <- if (filtered) "CEC: retained regions" else "CEC: complete diagram"
  graphics::plot(NA_real_, NA_real_, xlim = range(c(tiles$u0, tiles$u1)),
    ylim = range(c(tiles$lambda0, tiles$lambda1)), xlab = "log(C)", ylab = expression(lambda), main = main)
  graphics::rect(tiles$u0, tiles$lambda0, tiles$u1, tiles$lambda1, col = fill, border = NA)
  graphics::abline(h = 1, lty = 2)
  measured <- which(tiles$measured)
  graphics::points(log(tiles$C[measured]), tiles$lambda[measured], pch = 3, cex = .6)
  graphics::legend("topright", legend = c(paste(K, "groups"), if (filtered) "Excluded"),
    fill = c(colors, if (filtered) "#E1E4E8"), cex = .75, bg = "white")
  invisible(fit)
}

bands_validate <- function(fit) {
  if (!inherits(fit, "cec_bands")) stop("Supply a cec_bands result.")
  stopifnot(length(fit$runs) == fit$L + fit$B,
    all(vapply(fit$runs[seq_len(fit$L)], `[[`, "", "kind") == "original"),
    all(vapply(fit$runs[fit$L + seq_len(fit$B)], `[[`, "", "kind") == "bootstrap"),
    !anyDuplicated(vapply(fit$runs, `[[`, 0L, "seed")),
    all(fit$tiles$saturation[fit$tiles$keep] <= fit$control$saturation_max))
  for (run in fit$runs) {
    bands_check_indices(run$inner_indices, nrow(fit$X))
    for (candidate in run$candidates) bands_load_labels(candidate, fit$store)
  }
  for (i in seq_len(nrow(fit$points))) {
    result <- readRDS(file.path(fit$store$path, fit$points$checkpoint[i]))$value
    if (!result$valid) next
    allocations <- lapply(result$winners$allocation_key, bands_read_allocation, store = fit$store)
    recalculated <- bands_agreements(allocations, result$winners$allocation_key, fit$L, fit$B, bands_cache(0))
    stopifnot(identical(recalculated, result$diagnostics),
      result$diagnostics$original_pairs == choose(fit$L, 2),
      result$diagnostics$bootstrap_pairs == fit$L * fit$B, !result$conditional)
  }
  list(passed = TRUE, banks = length(fit$runs), points = nrow(fit$points),
    retained = nrow(bands_partitions(fit)), limits = fit$limits)
}

bands_read_allocation <- function(key, store) {
  labels <- bands_unpack_labels(readRDS(file.path(store$path, "allocations", paste0(key, ".rds"))))
  if (!identical(digest::digest(labels, algo = "xxhash64"), key)) stop("Allocation checksum mismatch.")
  labels
}
