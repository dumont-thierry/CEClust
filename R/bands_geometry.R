bands_entropy <- function(candidate) -sum(candidate$stats$mass * log(candidate$stats$mass))
bands_regular <- function(candidate) all(candidate$stats$rank == candidate$d)

bands_peaks <- function(candidate) {
  st <- candidate$stats
  exp(-candidate$d / 2 * log(2 * pi) - st$logdet[st$rank == candidate$d] / 2)
}

bands_lines <- function(candidates, C, penalty) {
  values <- lapply(candidates, bands_candidate_score, C = C, lambda = 0, penalty = penalty)
  list(A = vapply(values, `[[`, 0, "J"), H = vapply(values, `[[`, 0, "entropy"),
       saturation = vapply(values, `[[`, 0, "saturation"))
}

bands_candidate_score <- function(candidate, C, lambda, penalty) {
  bands_profile(candidate$stats, C, lambda, candidate$d, penalty)
}

bands_hull <- function(A, H, limits) {
  ordered <- order(-H, A, seq_along(A)); hull <- integer(); entry <- numeric()
  for (i in ordered) {
    if (length(hull) && H[i] == H[tail(hull, 1L)]) next
    crossing <- -Inf
    while (length(hull)) {
      previous <- tail(hull, 1L)
      crossing <- (A[i] - A[previous]) / (H[previous] - H[i])
      if (crossing > tail(entry, 1L)) break
      hull <- head(hull, -1L); entry <- head(entry, -1L)
    }
    if (!length(hull)) crossing <- -Inf
    hull <- c(hull, i); entry <- c(entry, crossing)
  }
  lower <- pmax(entry, limits[1L]); upper <- pmin(c(tail(entry, -1L), Inf), limits[2L])
  keep <- which(lower < upper)
  data.frame(partition = hull[keep], lower = lower[keep], upper = upper[keep])
}

bands_window <- function(candidates) {
  one <- which(vapply(candidates, `[[`, 0L, "K") == 1L)[1L]
  if (is.na(one)) stop("The reference pool requires a trivial partition.")
  C0 <- bands_peaks(candidates[[one]])
  peaks <- unlist(lapply(candidates, bands_peaks), use.names = FALSE)
  if (!length(peaks) || any(!is.finite(peaks) | peaks <= 0))
    stop("A density scale exceeds numerical range; no silent C truncation is applied.")
  C <- c(min(C0 / 8, min(peaks) / 2), max(C0 * 2^10, peaks * (1 + 1e-8)))
  if (any(!is.finite(C) | C <= 0)) stop("The calibrated C window exceeds numerical range.")
  list(C0 = C0, C = C, peaks = peaks, one = one)
}

bands_thin_centers <- function(values, maximum) {
  values <- sort(unique(values))
  if (length(values) <= maximum) return(values)
  values[unique(round(seq(1, length(values), length.out = maximum)))]
}

bands_make_tiles <- function(candidates, centers, limits, penalty) {
  u <- log(centers); mids <- (head(u, -1L) + tail(u, -1L)) / 2
  edges <- c(u[1L], mids, tail(u, 1L)); columns <- vector("list", length(centers))
  for (j in seq_along(centers)) {
    line <- bands_lines(candidates, centers[j], penalty)
    region <- bands_hull(line$A, line$H, limits); ids <- region$partition
    lambda <- (region$lower + region$upper) / 2
    columns[[j]] <- data.frame(column = j, u0 = edges[j], u1 = edges[j + 1L],
      C = centers[j], lambda = lambda, lambda0 = region$lower, lambda1 = region$upper,
      partition = ids, key = names(candidates)[ids],
      K = vapply(candidates[ids], `[[`, 0L, "K"),
      J = line$A[ids] + lambda * line$H[ids], saturation = line$saturation[ids])
  }
  result <- do.call(rbind, columns); rownames(result) <- NULL
  result$id <- seq_len(nrow(result)); result
}

bands_changed_columns <- function(tiles) {
  band <- ifelse(tiles$saturation == 0, 0L, ifelse(tiles$saturation <= .1, 1L, 2L))
  columns <- split(paste(tiles$key, band), tiles$column)
  changed <- logical(length(columns) - 1L)
  for (j in seq_along(changed)) changed[j] <- !identical(columns[[j]], columns[[j + 1L]])
  which(changed)
}

bands_geometry <- function(candidates, control, penalty) {
  started <- proc.time()[3L]; window <- bands_window(candidates)
  line <- bands_lines(candidates, window$C[2L], penalty)
  nontrivial <- line$H > 0
  threshold <- max(c(0, (line$A[window$one] - line$A[nontrivial]) / line$H[nontrivial]))
  limits <- c(0, max(1.5, threshold + .05 * max(1, abs(threshold))))
  releases <- vapply(candidates[vapply(candidates, bands_regular, TRUE)],
    bands_max_peak, 0)
  centers <- sort(unique(c(exp(seq(log(window$C[1L]), log(window$C[2L]),
    length.out = control$C_points)), releases * (1 - 1e-8), releases * (1 + 1e-8))))
  centers <- centers[centers >= window$C[1L] & centers <= window$C[2L]]
  search_limits <- character()
  if (length(centers) > control$max_C_points) search_limits <- "critical_C_sampling_budget"
  centers <- bands_thin_centers(centers, control$max_C_points)
  tiles <- bands_make_tiles(candidates, centers, limits, penalty); history <- list()
  for (round in seq_len(control$refine_rounds)) {
    changed <- bands_changed_columns(tiles)
    if (!length(changed)) break
    available <- control$max_C_points - length(centers)
    if (available <= 0L) { search_limits <- c(search_limits, "C_refinement_budget"); break }
    added <- sqrt(centers[changed] * centers[changed + 1L])
    extreme <- !is.finite(added) | added == 0
    if (any(extreme)) added[extreme] <- exp((log(centers[changed[extreme]]) +
      log(centers[changed[extreme] + 1L])) / 2)
    if (length(added) > available) {
      added <- bands_thin_centers(added, available)
      search_limits <- c(search_limits, "C_refinement_budget")
    }
    centers <- sort(unique(c(centers, added)))
    tiles <- bands_make_tiles(candidates, centers, limits, penalty)
    history[[round]] <- data.frame(round = round, added = length(added), cells = nrow(tiles))
  }
  if (length(bands_changed_columns(tiles))) search_limits <- c(search_limits, "finite_C_refinement")
  window$lambda <- limits; window$peaks <- NULL
  window$scope <- "Known candidates only. Left edge fully saturated; all known regular release scales lie before the right edge."
  list(tiles = tiles, C = centers, window = window, history = history,
    limits = unique(search_limits), seconds = unname(proc.time()[3L] - started))
}

bands_max_peak <- function(candidate) max(c(0, bands_peaks(candidate)))

bands_components <- function(tiles) {
  parent <- seq_len(nrow(tiles)); columns <- split(seq_len(nrow(tiles)), tiles$column)
  if (length(columns) > 1L) for (column in seq_len(length(columns) - 1L)) {
    a <- columns[[column]]; b <- columns[[column + 1L]]
    for (key in intersect(tiles$key[a], tiles$key[b])) {
      i <- a[match(key, tiles$key[a])]; j <- b[match(key, tiles$key[b])]
      overlap <- min(tiles$lambda1[i], tiles$lambda1[j]) - max(tiles$lambda0[i], tiles$lambda0[j])
      if (overlap > 0) parent[bands_find(j, parent)] <- bands_find(i, parent)
    }
  }
  roots <- vapply(seq_len(nrow(tiles)), bands_find, 0L, parent = parent)
  match(roots, unique(roots))
}

bands_find <- function(i, parent) {
  while (parent[i] != i) i <- parent[i]
  i
}

bands_order_component <- function(tiles, ids, saturation_max) {
  area <- (tiles$u1[ids] - tiles$u0[ids]) * (tiles$lambda1[ids] - tiles$lambda0[ids])
  u <- stats::weighted.mean(log(tiles$C[ids]), area)
  lambda <- stats::weighted.mean(tiles$lambda[ids], area)
  eligible <- ids[tiles$saturation[ids] <= saturation_max]
  distance <- (log(tiles$C[eligible]) - u)^2 + (tiles$lambda[eligible] - lambda)^2
  list(indices = eligible[order(distance, tiles$id[eligible])], logC = u, lambda = lambda)
}

bands_choose_anchors <- function(pilot, control) {
  tiles <- pilot$diagram$tiles; components <- bands_components(tiles)
  ids <- split(seq_len(nrow(tiles)), components); chosen <- integer()
  for (group in ids) {
    ordered <- bands_order_component(tiles, group, control$measurement_saturation)$indices
    ordered <- ordered[tiles$K[ordered] < pilot$n]
    if (length(ordered)) chosen <- c(chosen, ordered[1L])
  }
  # Fall back to a trivial anchor if no low-saturation component is available.
  if (!length(chosen)) chosen <- which(tiles$K == 1L)[1L]
  if (is.na(chosen[1L])) stop("No supported pilot anchor.")
  # Thin bands receive priority; the budget is a reported exploration limit.
  width <- tiles$lambda1[chosen] - tiles$lambda0[chosen]
  chosen <- chosen[order(width, tiles$K[chosen], tiles$id[chosen])]
  limited <- length(chosen) > control$max_anchors
  chosen <- head(chosen, control$max_anchors)
  result <- tiles[chosen, c("C", "lambda", "key", "K", "saturation", "lambda0", "lambda1")]
  result$id <- seq_len(nrow(result)); rownames(result) <- NULL
  attr(result, "limited") <- limited; result
}
