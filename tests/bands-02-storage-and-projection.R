library(CEClust)
N <- asNamespace("CEClust")

# Independent dense Gaussian projection used only for these small fixtures.
dense_projection <- function(training, labels, query, C, lambda) {
  labels <- N$bands_labels(labels)
  groups <- split(seq_len(nrow(training)), labels)
  values <- matrix(NA_real_, nrow(query), length(groups))
  for (k in seq_along(groups)) {
    sample <- training[groups[[k]], , drop = FALSE]
    center <- colMeans(sample)
    residual <- sweep(sample, 2, center)
    covariance <- crossprod(residual) / nrow(sample)
    determinant <- as.numeric(determinant(covariance, logarithm = TRUE)$modulus)
    bound <- -ncol(query) * log(2 * pi) - 2 * log(C)
    covariance <- covariance * exp(max(0, (bound - determinant) / ncol(query)))
    residual <- sweep(query, 2, center)
    distance <- rowSums((residual %*% solve(covariance)) * residual)
    values[, k] <- -.5 * (ncol(query) * log(2 * pi) +
      as.numeric(determinant(covariance, logarithm = TRUE)$modulus) + distance) +
      lambda * log(nrow(sample) / nrow(training))
  }
  N$bands_labels(max.col(values, ties.method = "first"))
}

directory <- tempfile("bands-storage-test-")
store <- N$bands_store(directory)

# Every encoding is lossless, including thresholds and partial final bytes.
encodings <- list(
  trivial = rep(1L, 17),
  singletons = seq_len(17),
  bits = rep(c(1L, 2L), length.out = 17),
  uint8 = rep(seq_len(255), 2),
  uint16le = rep(seq_len(256), 2),
  integer = rep(seq_len(65536), 2))
for (name in names(encodings)) {
  packed <- N$bands_pack_labels(encodings[[name]])
  stopifnot(identical(packed$encoding, name),
            identical(N$bands_unpack_labels(packed), encodings[[name]]))
}
stopifnot(length(N$bands_pack_labels(encodings$bits)$values) == 3L)

# Cache eviction changes availability, never the declared memory budget.
small_cache <- N$bands_cache(1024)
for (i in seq_len(20)) {
  N$bands_cache_put(small_cache, paste0("value-", i), rep(i, 50))
  stopifnot(small_cache$bytes <= small_cache$max_bytes)
}
stopifnot(small_cache$evictions > 0L)

X <- matrix(c(-.4, 0, .4, 2.6, 2.8, 3, 3.2, 3.4), ncol = 1)
labels <- c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 2L)
a <- N$bands_candidate(X, labels, store, "fixture")
b <- N$bands_candidate(X, ifelse(labels == 1L, 22L, 11L), store, "renumbered")
stopifnot(identical(a$key, b$key), is.null(a$labels),
          identical(N$bands_load_labels(a, store), labels))
stopifnot(identical(a$stats, N$bands_stats(X, labels)),
          identical(a$stats, b$stats), identical(a$rule, b$rule))
archive <- readRDS(file.path(store$path, a$label_file))
stopifnot(identical(archive$moments_format, 1L))

# The same labels fitted to another sample are not the same fitted rule.
c <- N$bands_candidate(X + 1, labels, store, "other-sample")
stopifnot(!identical(a$key, c$key))

# Both implicit endpoint encodings recover the exact hard partition.
one <- N$bands_candidate(X, rep(1L, nrow(X)), store, "trivial")
singletons <- N$bands_candidate(X, seq_len(nrow(X)), store, "endpoint")
stopifnot(identical(N$bands_load_labels(one, store), rep(1L, nrow(X))),
          identical(N$bands_load_labels(singletons, store), seq_len(nrow(X))))
endpoint_projection <- N$bands_project(singletons, X, 1, 1, store)
stopifnot(isTRUE(endpoint_projection$supported),
          identical(endpoint_projection$labels, seq_len(nrow(X))))
missing_point <- N$bands_candidate(X[-1, , drop = FALSE], seq_len(nrow(X) - 1L),
                                   store, "bootstrap-endpoint")
missing_projection <- N$bands_project(missing_point, X, 1, 1, store)
stopifnot(!isTRUE(missing_projection$supported),
          identical(missing_projection$unsupported, 1L))

query <- matrix(seq(-.6, 3.6, length.out = 31), ncol = 1)
cache <- N$bands_cache(1024^2)
for (ceiling in c(.1, 2, Inf)) for (lambda in c(.2, 1, 2)) {
  projected <- N$bands_project(a, query, ceiling, lambda, store,
                               cache = cache, target_key = "fixture-query")
  expected <- dense_projection(X, labels, query, ceiling, lambda)
  stopifnot(isTRUE(projected$supported), identical(projected$labels, expected))
}

# Reusing a certified lambda interval must agree with a fresh projection.
center <- N$bands_project(a, query, 2, 1, store, cache = cache,
                          target_key = "certificate-query")
interval <- center$interval
stopifnot(length(interval) == 2L, interval[1] <= 1, interval[2] >= 1)
inside <- (max(0, interval[1]) + min(2, interval[2])) / 2
if (inside > interval[1] && inside < interval[2] && inside > 0) {
  reused <- N$bands_project(a, query, 2, inside, store, cache = cache,
                            target_key = "certificate-query")
  fresh <- N$bands_project(a, query, 2, inside, store)
  stopifnot(identical(reused$labels, fresh$labels),
            identical(reused$labels, center$labels))
}

# A different target matrix must not hit a stale allocation cache entry.
changed_query <- query + .9
changed <- N$bands_project(a, changed_query, 2, 1, store, cache = cache)
stopifnot(identical(changed$labels,
                   dense_projection(X, labels, changed_query, 2, 1)))

# Singular supports are not extended by hidden regularization.
square <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
singular <- N$bands_candidate(square, c(1, 1, 2, 2), store, "singular")
unsupported <- N$bands_project(singular, matrix(c(.5, .5), nrow = 1), 1, 1,
                               store)
stopifnot(!isTRUE(unsupported$supported), is.null(unsupported$labels),
          identical(as.integer(unsupported$unsupported), 1L))

# A fresh store handle must retrieve candidates archived by an earlier handle.
reopened <- N$bands_store(directory)
stopifnot(identical(N$bands_load_labels(a, reopened), labels))
