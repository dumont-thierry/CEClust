# Gaussian numerical definitions adapted from CEClustNext (Thierry Dumont, MIT).
# Finite C uses the same raw determinant-constrained criterion and CEM proposal.

bands_data <- function(X) {
  if (is.numeric(X) && is.null(dim(X))) X <- matrix(X, ncol = 1L)
  X <- as.matrix(X)
  if (!is.numeric(X) || nrow(X) < 2L || ncol(X) < 1L || any(!is.finite(X)))
    stop("X must be a finite numeric matrix with at least two rows.")
  storage.mode(X) <- "double"
  X
}

bands_labels <- function(labels) {
  if (!is.atomic(labels) || !is.null(dim(labels)) || !length(labels) || anyNA(labels))
    stop("Labels must be a nonempty vector without missing values.")
  as.integer(match(labels, unique(labels)))
}

bands_seed <- function(...) {
  strtoi(substr(digest::digest(list(...), algo = "xxhash64"), 1L, 7L), 16L)
}

bands_numeric_integer <- function(value, name, minimum = 1L, maximum = .Machine$integer.max) {
  if (length(value) != 1L || !is.numeric(value) || !is.finite(value) ||
      value != floor(value) || value < minimum || value > maximum)
    stop(name, " must be an integer between ", minimum, " and ", maximum, ".")
  as.integer(value)
}

bands_initial <- function(X, order, shot) {
  X <- bands_data(X); n <- nrow(X)
  order <- bands_numeric_integer(order, "order", 1L, n)
  shot <- bands_numeric_integer(shot, "shot")
  if (order == 1L) return(rep(1L, n))
  if (shot %% 5L == 0L) return(sample(rep(seq_len(order), length.out = n)))
  centers <- sample.int(n, 1L); distance <- rep(Inf, n)
  for (j in seq.int(2L, order)) {
    difference <- sweep(X, 2L, X[utils::tail(centers, 1L), ], "-")
    distance <- pmin(distance, rowSums(difference * difference)); distance[centers] <- 0
    remaining <- setdiff(seq_len(n), centers)
    next_id <- if (sum(distance) > 0) sample.int(n, 1L, prob = distance) else
      remaining[sample.int(length(remaining), 1L)]
    centers <- c(centers, next_id)
  }
  fit <- tryCatch(suppressWarnings(stats::kmeans(X, centers = X[centers, , drop = FALSE],
      iter.max = 50L, algorithm = "Lloyd")), error = identity)
  labels <- if (inherits(fit, "error")) sample(rep(seq_len(order), length.out = n)) else fit$cluster
  if (length(unique(labels)) != order) labels <- sample(rep(seq_len(order), length.out = n))
  bands_labels(labels)
}

bands_stats <- function(X, labels) {
  X <- bands_data(X); labels <- bands_labels(labels)
  if (length(labels) != nrow(X)) stop("One label is required per observation.")
  groups <- split(seq_len(nrow(X)), labels); d <- ncol(X)
  result <- vector("list", length(groups))
  for (k in seq_along(groups)) {
    A <- X[groups[[k]], , drop = FALSE]
    shifted <- sweep(A, 2L, A[1L, ], "-")
    centered <- sweep(shifted, 2L, colMeans(shifted), "-")
    covariance <- crossprod(centered) / nrow(A)
    values <- eigen(covariance, symmetric = TRUE, only.values = TRUE)$values
    rank <- sum(values > 8 * d * .Machine$double.eps * max(abs(values)))
    determinant <- if (rank == d) 2 * sum(log(diag(chol(covariance)))) else -Inf
    result[[k]] <- data.frame(cluster = k, n = nrow(A), mass = nrow(A) / nrow(X),
                             rank = rank, logdet = determinant)
  }
  do.call(rbind, result)
}

bands_profile <- function(stats, C, lambda, d, penalty = "none") {
  penalty <- match.arg(penalty, c("none", "bic"))
  if (length(C) != 1L || !is.numeric(C) || is.na(C) || C <= 0 ||
      length(lambda) != 1L || !is.numeric(lambda) || !is.finite(lambda) || lambda < 0)
    stop("C must be positive (possibly Inf), and lambda finite and nonnegative.")
  d <- bands_numeric_integer(d, "d")
  n <- sum(stats$n); singular <- stats$rank < d
  if (!is.finite(C) && any(singular)) stop("An unbounded Gaussian profile requires full-rank groups.")
  if (is.finite(C)) {
    bound <- -d * log(2 * pi) - 2 * log(C)
    scale <- pmax(0, (bound - stats$logdet) / d)
    conditional <- d / 2 * log(2 * pi) + pmax(stats$logdet, bound) / 2 + d / 2 * exp(-scale)
    conditional[singular] <- -log(C)
    saturated <- singular | stats$logdet <= bound + 1e-10 * max(1, abs(bound))
  } else {
    conditional <- d / 2 * (log(2 * pi) + 1) + stats$logdet / 2
    saturated <- rep(FALSE, nrow(stats))
  }
  entropy <- -sum(stats$mass * log(stats$mass))
  beta <- if (penalty == "bic") (d + 1) * (d + 2) * log(n) / (4 * n) else 0
  list(J = lambda * entropy + sum(stats$mass * conditional) + beta * (nrow(stats) - 1),
       entropy = entropy, K = nrow(stats), saturation = sum(stats$mass[saturated]),
       penalty = beta * (nrow(stats) - 1))
}

bands_C0 <- function(X) {
  X <- bands_data(X); one <- bands_stats(X, rep(1L, nrow(X)))
  if (one$rank < ncol(X)) stop("The complete sample must have full covariance rank.")
  value <- exp(-ncol(X) / 2 * log(2 * pi) - one$logdet / 2)
  if (!is.finite(value) || value <= 0) stop("The trivial density peak is outside numerical range.")
  value
}

bands_fit <- function(X, labels, C, lambda, penalty, control) {
  X <- bands_data(X); labels <- bands_labels(labels)
  if (length(labels) != nrow(X)) stop("One label is required per observation.")
  penalty <- match.arg(penalty, c("none", "bic"))
  if (length(C) != 1L || !is.numeric(C) || is.na(C) || C <= 0 ||
      length(lambda) != 1L || !is.numeric(lambda) || !is.finite(lambda) || lambda <= 0)
    stop("C must be positive (possibly Inf), and lambda finite and positive.")
  iterations <- bands_numeric_integer(control$max_iter, "control$max_iter")
  tolerance <- control$tolerance
  if (length(tolerance) != 1L || !is.numeric(tolerance) || !is.finite(tolerance) || tolerance < 0)
    stop("control$tolerance must be finite and nonnegative.")
  d <- ncol(X); n <- nrow(X)
  beta <- if (penalty == "bic") (d + 1) * (d + 2) * log(n) / (4 * n) else 0
  fit <- tryCatch(bands_cem_cpp(X, labels, C, lambda, iterations, tolerance, beta), error = identity)
  if (inherits(fit, "error")) return(list(labels = labels, J = NA_real_, trace = numeric(),
    stop = paste0("numerical_error: ", conditionMessage(fit)), valid = FALSE))
  if (isTRUE(fit$valid) && (any(!is.finite(fit$trace)) || any(diff(fit$trace) >= 0)))
    stop("Accepted CEM steps must strictly decrease the finite criterion.")
  fit$iterations <- max(0L, length(fit$trace) - 1L)
  fit
}
