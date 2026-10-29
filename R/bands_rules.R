# Compact Gaussian allocation rules. The density/support convention follows
# CEClustNext 0.1.0.9000; block evaluation does not retain an n-by-K matrix.
bands_rule_group <- function(X) {
  d <- ncol(X); shifted <- sweep(X, 2L, X[1L, ], "-")
  offset <- colMeans(shifted); centered <- sweep(shifted, 2L, offset, "-")
  covariance <- crossprod(centered) / nrow(X)
  values <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
  tolerance <- 8 * d * .Machine$double.eps * max(abs(values))
  rank <- sum(values > tolerance)
  list(mean = X[1L, ] + offset, covariance = covariance, rank = rank,
       logdet = if (rank == d) 2 * sum(log(diag(chol(covariance)))) else -Inf)
}

bands_rule <- function(X, labels) {
  X <- bands_data(X); labels <- bands_labels(labels)
  if (length(labels) != nrow(X)) stop("Allocation and training sample differ in length.")
  if (max(labels) == nrow(X)) return(list(n = nrow(X), d = ncol(X), K = nrow(X),
    type = "singletons", points = X, regular = FALSE, ceiling = Inf))
  indices <- split(seq_len(nrow(X)), labels)
  groups <- lapply(indices, bands_rule_rows, data = X)
  count <- lengths(indices); logmass <- log(count / nrow(X))
  regular <- vapply(groups, bands_group_regular, TRUE, d = ncol(X))
  peaks <- vapply(groups[regular], bands_group_peak, 0, d = ncol(X))
  list(n = nrow(X), d = ncol(X), K = length(groups), type = "gaussian", groups = groups,
       count = unname(count), logmass = logmass, regular = all(regular),
       ceiling = if (all(regular)) max(peaks) else Inf)
}

bands_rule_stats <- function(rule) {
  if (identical(rule$type, "singletons")) return(data.frame(cluster = seq_len(rule$n),
    n = rep(1L, rule$n), mass = rep(1 / rule$n, rule$n),
    rank = integer(rule$n), logdet = rep(-Inf, rule$n)))
  data.frame(cluster = seq_len(rule$K), n = rule$count, mass = rule$count / rule$n,
    rank = unname(vapply(rule$groups, `[[`, 0L, "rank")),
    logdet = unname(vapply(rule$groups, `[[`, 0, "logdet")))
}

bands_rule_rows <- function(rows, data) bands_rule_group(data[rows, , drop = FALSE])
bands_group_regular <- function(group, d) group$rank == d
bands_group_peak <- function(group, d) exp(-d / 2 * log(2 * pi) - group$logdet / 2)

bands_prepare_group <- function(group, C, d) {
  covariance <- group$covariance
  if (group$rank == d) {
    bound <- -d * log(2 * pi) - 2 * log(C)
    covariance <- covariance * exp(max(0, (bound - group$logdet) / d))
  }
  decomposition <- eigen(covariance, symmetric = TRUE)
  tolerance <- 8 * d * .Machine$double.eps * max(abs(decomposition$values))
  active <- decomposition$values > tolerance
  if (all(active)) {
    L <- chol(covariance)
    return(list(mean = group$mean, regular = TRUE, L = L,
                constant = -d / 2 * log(2 * pi) - sum(log(diag(L)))))
  }
  if (!is.finite(C)) stop("An unbounded allocation rule cannot contain singular groups.")
  list(mean = group$mean, regular = FALSE,
       null_basis = decomposition$vectors[, !active, drop = FALSE],
       tolerance = tolerance, logC = log(C))
}

bands_group_density <- function(group, X) {
  centered <- sweep(X, 2L, group$mean, "-")
  if (group$regular) {
    solved <- backsolve(group$L, t(centered), transpose = TRUE)
    return(group$constant - colSums(solved * solved) / 2)
  }
  residual <- sqrt(rowSums((centered %*% group$null_basis)^2))
  allowance <- sqrt(max(0, group$tolerance)) + 64 * ncol(X) * .Machine$double.eps *
    pmax(apply(abs(X), 1L, max), max(abs(group$mean)))
  ifelse(residual <= allowance, group$logC, -Inf)
}

bands_rule_unbounded <- function(rule, C) {
  if (!rule$regular) return(FALSE)
  bound <- -rule$d * log(2 * pi) - 2 * log(C)
  all(vapply(rule$groups, bands_group_released, TRUE, bound = bound))
}

bands_group_released <- function(group, bound) bound <= group$logdet

bands_projection_interval <- function(density, choices, logmass, limits) {
  rows <- seq_len(nrow(density)); selected <- density[cbind(rows, choices)]
  chosen_mass <- logmass[choices]
  for (k in seq_along(logmass)) {
    usable <- is.finite(density[, k])
    difference <- chosen_mass[usable] - logmass[k]
    gap <- selected[usable] - density[usable, k]
    positive <- difference > 0; negative <- difference < 0
    if (any(positive)) limits[1L] <- max(limits[1L], max(-gap[positive] / difference[positive]))
    if (any(negative)) limits[2L] <- min(limits[2L], min(-gap[negative] / difference[negative]))
  }
  limits
}

bands_cached_projection <- function(entries, lambda) {
  if (is.null(entries)) return(NULL)
  for (entry in entries) if (isTRUE(lambda == entry$lambda) ||
      (lambda > entry$interval[1L] && lambda < entry$interval[2L])) return(entry)
  NULL
}

bands_singleton_projection <- function(rule, Xnew, C, lambda) {
  if (!is.finite(C)) stop("An unbounded allocation rule cannot contain singular groups.")
  means <- rule$points
  distinct <- which(!duplicated(means))
  ordered <- distinct[order(means[distinct, 1L])]
  first <- means[ordered, 1L]
  allowance <- 64 * rule$d * .Machine$double.eps * max(abs(Xnew), abs(means))
  if (length(first) > 1L && any(diff(first) <= 2 * allowance)) return(NULL)
  interval <- findInterval(Xnew[, 1L], first)
  left <- pmax(1L, pmin(length(first), interval))
  right <- pmax(1L, pmin(length(first), interval + 1L))
  use_right <- abs(Xnew[, 1L] - first[right]) < abs(Xnew[, 1L] - first[left])
  selected <- ordered[ifelse(use_right, right, left)]
  differences <- Xnew - means[selected, , drop = FALSE]
  tolerance <- 64 * rule$d * .Machine$double.eps *
    pmax(apply(abs(Xnew), 1L, max), apply(abs(means[selected, , drop = FALSE]), 1L, max))
  unsupported <- which(sqrt(rowSums(differences * differences)) > tolerance)
  list(labels = if (length(unsupported)) NULL else bands_labels(selected),
       supported = length(unsupported) == 0L, unsupported = as.integer(unsupported),
       interval = c(0, Inf), lambda = lambda, route = "singleton_coordinate_lookup")
}

bands_singleton_density <- function(rule, X, C, group) {
  centered <- sweep(X, 2L, rule$points[group, ], "-")
  tolerance <- 64 * rule$d * .Machine$double.eps *
    pmax(apply(abs(X), 1L, max), max(abs(rule$points[group, ])))
  ifelse(sqrt(rowSums(centered * centered)) <= tolerance, log(C), -Inf)
}

bands_guard_interval <- function(limits, lambda) {
  finite <- is.finite(limits)
  margin <- 64 * .Machine$double.eps * pmax(1, abs(limits[finite]))
  limits[finite] <- limits[finite] + c(1, -1)[finite] * margin
  if (lambda <= limits[1L] || lambda >= limits[2L]) c(lambda, lambda) else limits
}

bands_projection_data <- function(X) {
  if (is.numeric(X) && is.null(dim(X))) X <- matrix(X, ncol = 1L)
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) < 1L || ncol(X) < 1L || any(!is.finite(X)))
    stop("Projection requires a nonempty finite numeric matrix.")
  storage.mode(X) <- "double"
  X
}

bands_project <- function(candidate, Xnew, C, lambda, store = NULL, cache = NULL,
                          target_key = NULL, block_bytes = 8 * 1024^2) {
  Xnew <- bands_projection_data(Xnew); rule <- candidate$rule
  if (ncol(Xnew) != rule$d) stop("Projection dimension mismatch.")
  if (length(C) != 1L || !is.numeric(C) || is.na(C) || C <= 0 ||
      length(lambda) != 1L || !is.numeric(lambda) || !is.finite(lambda) || lambda <= 0)
    stop("Projection requires positive C and finite positive lambda.")
  if (length(block_bytes) != 1L || !is.finite(block_bytes) || block_bytes <= 0)
    stop("Positive block budget required.")
  masses <- if (identical(rule$type, "singletons")) -log(rule$n) else rule$logmass
  if (any(!is.finite(lambda * masses))) stop("Projection allocation weights overflow.")
  if (is.null(target_key)) target_key <- digest::digest(Xnew, algo = "sha256")
  ceiling_key <- if (bands_rule_unbounded(rule, C)) "unbounded" else sprintf("%.17g", C)
  key <- paste("projection", candidate$key, target_key, ceiling_key, sep = ":")
  entries <- bands_cache_get(cache, key); old <- bands_cached_projection(entries, lambda)
  if (!is.null(old)) {
    old$evaluated_lambda <- old$lambda; old$lambda <- lambda
    old$route <- "certified_interval_cache"; return(old)
  }
  singleton <- identical(rule$type, "singletons")
  if (singleton) {
    value <- bands_singleton_projection(rule, Xnew, C, lambda)
    if (!is.null(value)) {
      entries[[length(entries) + 1L]] <- value
      bands_cache_put(cache, key, entries)
      return(value)
    }
  }
  prepared <- if (singleton) NULL else lapply(rule$groups, bands_prepare_group, C = C, d = rule$d)
  logmass <- if (singleton) rep(-log(rule$n), rule$K) else rule$logmass
  block_rows <- max(1L, min(nrow(Xnew), floor(block_bytes / (24 * rule$K))))
  choices <- integer(nrow(Xnew)); unsupported <- integer(); limits <- c(0, Inf)
  for (first in seq.int(1L, nrow(Xnew), by = block_rows)) {
    rows <- seq.int(first, min(nrow(Xnew), first + block_rows - 1L))
    block <- Xnew[rows, , drop = FALSE]
    density <- matrix(0, length(rows), rule$K)
    for (k in seq_len(rule$K)) density[, k] <- if (singleton)
      bands_singleton_density(rule, block, C, k) else
        bands_group_density(prepared[[k]], block)
    scores <- sweep(density, 2L, lambda * logmass, "+")
    supported <- rowSums(is.finite(scores)) > 0L
    unsupported <- c(unsupported, rows[!supported])
    if (!any(supported)) next
    local <- max.col(scores[supported, , drop = FALSE], ties.method = "first")
    choices[rows[supported]] <- local
    limits <- bands_projection_interval(density[supported, , drop = FALSE], local, logmass, limits)
  }
  limits <- if (length(unsupported)) c(lambda, lambda) else bands_guard_interval(limits, lambda)
  value <- list(labels = if (length(unsupported)) NULL else bands_labels(choices),
                supported = length(unsupported) == 0L, unsupported = as.integer(unsupported),
                interval = limits, lambda = lambda, route = "compact_rule_blocks")
  entries[[length(entries) + 1L]] <- value
  bands_cache_put(cache, key, entries)
  value
}
