bands_bank_winner <- function(bank, C, lambda, penalty, cache = NULL) {
  key <- paste("lines", bank$signature, format(C, digits = 17), penalty, sep = ":")
  line <- bands_cache_get(cache, key)
  if (is.null(line)) {
    line <- bands_lines(bank$candidates, C, penalty)
    bands_cache_put(cache, key, line)
  }
  score <- line$A + lambda * line$H; index <- which.min(score)
  list(candidate = bank$candidates[[index]], J = score[index], saturation = line$saturation[index])
}

bands_pool_originals <- function(runs) {
  candidates <- list(); provenance <- list()
  for (run in runs) if (run$kind == "original") for (candidate in run$candidates) {
    candidates <- bands_add_candidate(candidates, candidate)
    provenance[[length(provenance) + 1L]] <- data.frame(key = candidate$key,
      run = run$name, source_key = candidate$key, C = NA_real_, lambda = NA_real_,
      origin = "original_bank")
  }
  list(candidates = candidates, provenance = do.call(rbind, provenance))
}

bands_spread_indices <- function(n) {
  result <- integer(); pending <- list(c(1L, n)); cursor <- 1L
  while (cursor <= length(pending)) {
    interval <- pending[[cursor]]; cursor <- cursor + 1L
    middle <- as.integer(floor(mean(interval))); result <- c(result, middle)
    if (interval[1L] < middle) pending[[length(pending) + 1L]] <- c(interval[1L], middle - 1L)
    if (middle < interval[2L]) pending[[length(pending) + 1L]] <- c(middle + 1L, interval[2L])
  }
  result
}

bands_projection_queries <- function(anchors, diagram) {
  tiles <- diagram$tiles; groups <- split(seq_len(nrow(tiles)), bands_components(tiles))
  representatives <- integer()
  for (group in groups) {
    ordered <- bands_order_component(tiles, group, 1)$indices
    if (length(ordered)) representatives <- c(representatives, ordered[1L])
  }
  representatives <- representatives[order(tiles$saturation[representatives] > .1,
    tiles$saturation[representatives], tiles$K[representatives])]
  rest <- bands_spread_indices(nrow(tiles))
  rest <- rest[order(tiles$saturation[rest] > .1, seq_along(rest))]
  # All regions remain eligible; low-saturation component centers receive priority.
  points <- rbind(anchors[, c("C", "lambda")],
    tiles[representatives, c("C", "lambda")], tiles[rest, c("C", "lambda")])
  points <- unique(points); rownames(points) <- NULL; points
}

bands_reference <- function(X, runs, anchors, control, penalty, store) {
  started <- proc.time()[3L]; pooled <- bands_pool_originals(runs)
  candidates <- pooled$candidates; provenance <- list(pooled$provenance)
  data_key <- digest::digest(X); cache <- bands_cache(control$cache_mb * 1024^2)
  seen <- new.env(hash = TRUE, parent = emptyenv()); queries <- 0L; unsupported <- 0L
  history <- list(); limits <- character(); stable <- FALSE
  boots <- runs[vapply(runs, bands_is_bootstrap, TRUE)]
  diagram <- bands_geometry(candidates, control, penalty)
  for (round in seq_len(control$projection_rounds)) {
    points <- bands_projection_queries(anchors, diagram); before <- length(candidates)
    exhausted <- FALSE
    for (i in seq_len(nrow(points))) {
      q <- points[i, ]
      keys <- vapply(boots, bands_query_key, "", C = q$C, lambda = q$lambda)
      pending <- which(!vapply(keys, exists, TRUE, envir = seen, inherits = FALSE))
      if (queries + length(pending) > control$max_projection_queries) { exhausted <- TRUE; break }
      for (b in pending) {
        run <- boots[[b]]; key <- keys[b]
        assign(key, TRUE, seen); queries <- queries + 1L
        value <- bands_reference_projection(run, X, q$C, q$lambda, penalty, store, cache, data_key)
        if (is.null(value$candidate)) { unsupported <- unsupported + 1L; next }
        candidate <- value$candidate
        candidates <- bands_add_candidate(candidates, candidate)
        provenance[[length(provenance) + 1L]] <- data.frame(key = candidate$key,
          run = run$name, source_key = value$source_key, C = q$C, lambda = q$lambda,
          origin = "bootstrap_projection")
      }
    }
    diagram <- bands_geometry(candidates, control, penalty)
    history[[round]] <- data.frame(round = round, before = before, after = length(candidates),
      evaluated_projections = queries, unsupported = unsupported, cells = nrow(diagram$tiles))
    if (exhausted) { limits <- "projection_query_budget"; break }
    if (length(candidates) == before) { stable <- TRUE; break }
  }
  if (!stable) limits <- c(limits, "finite_projection_enrichment")
  list(candidates = candidates, diagram = diagram, provenance = do.call(rbind, provenance),
    history = history, queries = queries, unsupported = unsupported, stable = stable,
    limits = unique(c(limits, diagram$limits)), seconds = unname(proc.time()[3L] - started),
    scope = "Union of original banks and bootstrap winners projected at recorded queries. Not all possible projections on a continuum.")
}

bands_is_bootstrap <- function(run) identical(run$kind, "bootstrap")

bands_query_key <- function(run, C, lambda) {
  paste(run$name, sprintf("%.17g", C), sprintf("%.17g", lambda))
}

bands_reference_projection <- function(run, X, C, lambda, penalty, store, cache, data_key) {
  key <- digest::digest(list(run$signature, data_key, C, lambda, penalty))
  path <- file.path(store$path, "reference_queries", paste0(key, ".rds"))
  if (file.exists(path)) return(readRDS(path))
  winner <- bands_bank_winner(run, C, lambda, penalty, cache)
  projected <- bands_project(winner$candidate, X, C, lambda, store, cache, data_key)
  candidate <- if (is.null(projected$labels)) NULL else
    bands_candidate(X, projected$labels, store, "bootstrap_projection", data_key)
  result <- list(candidate = candidate, source_key = winner$candidate$key,
    C = C, lambda = lambda, supported = !is.null(candidate))
  bands_save(result, path); result
}
