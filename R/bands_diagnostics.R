bands_archive_allocation <- function(labels, store) {
  labels <- bands_labels(labels); key <- digest::digest(labels, algo = "xxhash64")
  path <- file.path(store$path, "allocations", paste0(key, ".rds"))
  if (file.exists(path)) {
    if (!identical(bands_unpack_labels(readRDS(path)), labels)) stop("Allocation hash collision.")
  } else bands_save(bands_pack_labels(labels), path)
  key
}

bands_matching <- function(a, b) {
  if (length(a) != length(b)) stop("Agreements require the same observations.")
  if (identical(a, b)) return(1)
  n <- length(a); Ka <- max(a); Kb <- max(b)
  if (Ka == n) return(Kb / n)
  if (Kb == n) return(Ka / n)
  if (Ka == 1L || Kb == 1L) return(max(tabulate(if (Ka == 1L) b else a)) / n)
  tab <- table(a, b); K <- max(dim(tab)); square <- matrix(0, K, K)
  square[seq_len(nrow(tab)), seq_len(ncol(tab))] <- tab
  assignment <- as.integer(clue::solve_LSAP(square, maximum = TRUE))
  sum(square[cbind(seq_len(K), assignment)]) / n
}

bands_pair_agreement <- function(a, b, ka, kb, cache) {
  if (identical(ka, kb)) return(1)
  key <- paste(c("agreement", sort(c(ka, kb))), collapse = ":")
  value <- bands_cache_get(cache, key)
  if (is.null(value)) { value <- bands_matching(a, b); bands_cache_put(cache, key, value) }
  value
}

bands_agreements <- function(labels, keys, L, B, cache) {
  key <- paste0("pattern:", digest::digest(list(keys, L, B)))
  previous <- bands_cache_get(cache, key)
  if (!is.null(previous)) return(previous)
  original <- numeric(choose(L, 2)); cross <- numeric(L * B); at <- 0L
  for (i in seq_len(L - 1L)) for (j in seq.int(i + 1L, L)) {
    at <- at + 1L
    original[at] <- bands_pair_agreement(labels[[i]], labels[[j]], keys[i], keys[j], cache)
  }
  at <- 0L
  for (i in seq_len(L)) for (b in seq_len(B)) {
    at <- at + 1L; j <- L + b
    cross[at] <- bands_pair_agreement(labels[[i]], labels[[j]], keys[i], keys[j], cache)
  }
  algo <- mean(original); boot <- mean(cross)
  result <- list(stab_algo = algo, stab_boot = boot, RSI = if (algo > 0) boot / algo else NA_real_,
    original_pairs = length(original), bootstrap_pairs = length(cross),
    original_agreements = original, bootstrap_agreements = cross)
  bands_cache_put(cache, key, result); result
}

bands_evaluate_point <- function(X, runs, candidates, C, lambda, L, B, control,
    penalty, store, cache) {
  data_key <- digest::digest(X); winners <- list(); labels <- list(); keys <- character()
  for (i in seq_along(runs)) {
    run <- runs[[i]]; winner <- bands_bank_winner(run, C, lambda, penalty, cache)
    if (run$kind == "original") {
      cache_key <- paste0("hard:", winner$candidate$key)
      allocation <- bands_cache_get(cache, cache_key)
      if (is.null(allocation)) {
        allocation <- bands_load_labels(winner$candidate, store)
        bands_cache_put(cache, cache_key, allocation)
      }
    } else allocation <- bands_project(winner$candidate, X, C, lambda, store, cache, data_key)$labels
    if (is.null(allocation)) return(list(valid = FALSE, passed = FALSE,
      reason = "unsupported_bootstrap_projection", failed_run = run$name, C = C, lambda = lambda))
    labels[[i]] <- allocation; keys[i] <- bands_archive_allocation(allocation, store)
    winners[[i]] <- data.frame(run = run$name, kind = run$kind,
      candidate_key = winner$candidate$key, allocation_key = keys[i], J = winner$J,
      K_fit = winner$candidate$K, K_projected = max(allocation))
  }
  diagnostics <- bands_agreements(labels, keys, L, B, cache)
  reference <- bands_lines(candidates, C, penalty)
  values <- reference$A + lambda * reference$H; best <- which.min(values)
  passed <- is.finite(diagnostics$RSI) && diagnostics$stab_algo > control$algo_min &&
    diagnostics$RSI > control$rsi_min
  list(valid = TRUE, passed = passed, C = C, lambda = lambda,
    reference_key = candidates[[best]]$key, reference_J = values[best],
    reference_saturation = reference$saturation[best], diagnostics = diagnostics,
    winners = do.call(rbind, winners), conditional = FALSE,
    scope = "Independent banks conditional on the shared pilot. No champion agreement gate or conditional RSI correction.")
}

bands_diagnose <- function(X, runs, reference, L, B, control, penalty, store, signature) {
  started <- proc.time()[3L]; tiles <- reference$diagram$tiles
  tiles$component <- bands_components(tiles)
  tiles$measured <- FALSE; tiles$witness <- NA_integer_; tiles$keep <- FALSE
  tiles$stab_algo <- tiles$stab_boot <- tiles$RSI <- NA_real_
  points <- list(); groups <- split(seq_len(nrow(tiles)), tiles$component)
  cache <- bands_cache(control$cache_mb * 1024^2)
  for (ids in groups) {
    ordered <- bands_order_component(tiles, ids, control$measurement_saturation)$indices
    for (i in ordered) {
      point <- tiles[i, ]; path <- file.path(store$path, "points", paste0(point$id, ".rds"))
      result <- NULL
      if (file.exists(path)) {
        saved <- readRDS(path)
        if (!identical(saved$signature, signature)) stop("Incompatible diagnostic checkpoint.")
        result <- saved$value
      }
      if (is.null(result)) {
        result <- bands_evaluate_point(X, runs, reference$candidates, point$C, point$lambda,
          L, B, control, penalty, store, cache)
        bands_save(list(signature = signature, value = result), path)
      }
      tiles$measured[i] <- TRUE
      points[[as.character(point$id)]] <- data.frame(id = point$id, C = point$C,
        lambda = point$lambda, valid = result$valid, passed = result$passed,
        stab_algo = if (result$valid) result$diagnostics$stab_algo else NA_real_,
        stab_boot = if (result$valid) result$diagnostics$stab_boot else NA_real_,
        RSI = if (result$valid) result$diagnostics$RSI else NA_real_,
        checkpoint = file.path("points", paste0(point$id, ".rds")))
      if (result$valid) {
        tiles$stab_algo[i] <- result$diagnostics$stab_algo
        tiles$stab_boot[i] <- result$diagnostics$stab_boot; tiles$RSI[i] <- result$diagnostics$RSI
      }
      if (!result$passed) next
      target <- if (control$diag_mode == "components") ids else i
      tiles$witness[target] <- point$id
      tiles$keep[target] <- tiles$saturation[target] <= control$saturation_max
      # Unmeasured cells inherit a decision, never fabricated numerical diagnostics.
      if (control$diag_mode == "components") break
    }
  }
  list(tiles = tiles, points = if (length(points)) do.call(rbind, points) else data.frame(),
    seconds = unname(proc.time()[3L] - started),
    cache = list(bytes = cache$bytes, hits = cache$hits, evictions = cache$evictions))
}
