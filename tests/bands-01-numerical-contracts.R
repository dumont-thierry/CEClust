library(CEClust)
N <- asNamespace("CEClust")

near <- function(actual, expected, tolerance = 1e-10) {
  stopifnot(isTRUE(all.equal(actual, expected, tolerance = tolerance,
                             check.attributes = FALSE)))
}

# The one-dimensional Gaussian profile has an independent closed form.
X <- matrix(c(-2, -1, 0, 1, 2), ncol = 1)
labels <- rep(1L, nrow(X))
variance <- mean((X[, 1] - mean(X[, 1]))^2)
st <- N$bands_stats(X, labels)
expected <- (log(2 * pi * variance) + 1) / 2
near(N$bands_profile(st, 10, 1, 1)$J, expected)
near(N$bands_C0(X), 1 / sqrt(2 * pi * variance))

# At a low density ceiling, the constrained covariance is 1/(2*pi*C^2).
C <- .01
saturated <- N$bands_profile(st, C, 1, 1)
near(saturated$J, -log(C) + pi * variance * C^2)
stopifnot(saturated$saturation == 1)

# A singleton endpoint is finite only with a finite density ceiling.
st_singletons <- N$bands_stats(X, seq_len(nrow(X)))
endpoint <- N$bands_profile(st_singletons, 3, 1.3, 1)
near(endpoint$J, -log(3) + 1.3 * log(nrow(X)))
stopifnot(endpoint$saturation == 1, endpoint$K == nrow(X))
unbounded_singular <- try(N$bands_profile(st_singletons, Inf, 1, 1), silent = TRUE)
stopifnot(inherits(unbounded_singular, "try-error"))

# Label names cannot change a hard-allocation score or its saturation.
X2 <- rbind(c(-2.1, -.1), c(-1.9, -.2), c(-2, .3),
            c(2.1, .2), c(1.9, .1), c(2, -.3))
a <- N$bands_stats(X2, c(1, 1, 1, 2, 2, 2))
b <- N$bands_stats(X2, c(99, 99, 99, 4, 4, 4))
near(N$bands_profile(a, 10, 1.2, 2),
     N$bands_profile(b, 10, 1.2, 2))
plain <- N$bands_profile(a, 10, 1.2, 2, "none")
penalized <- N$bands_profile(a, 10, 1.2, 2, "bic")
near(penalized$J - plain$J, 3 * log(6) / 6)

# A free descent must reject singular initial groups explicitly. The same
# starting allocations remain a defined finite-C problem.
square <- rbind(c(0, 0), c(1, 0), c(0, 1), c(1, 1))
control <- list(max_iter = 20L, tolerance = 1e-8)
free_bad <- N$bands_fit(square, c(1, 1, 2, 2), Inf, 1, "none", control)
bounded <- N$bands_fit(square, c(1, 1, 2, 2), 1, 1, "none", control)
stopifnot(!isTRUE(free_bad$valid), isTRUE(bounded$valid))
stopifnot(length(free_bad$stop) == 1L, nzchar(free_bad$stop))
stopifnot(all(is.finite(bounded$trace)),
          length(bounded$trace) < 2L || all(diff(bounded$trace) < 0))

# An ordinary regular descent reports its actual optimized score. This does
# not assume that it finds a global optimum or any prescribed group count.
initial <- c(1, 1, 1, 2, 2, 2)
for (ceiling in c(2, Inf)) {
  fit <- N$bands_fit(X2, initial, ceiling, 1.1, "none", control)
  stopifnot(isTRUE(fit$valid), length(fit$labels) == nrow(X2),
            all(is.finite(fit$trace)),
            length(fit$trace) < 2L || all(diff(fit$trace) < 0))
  final_stats <- N$bands_stats(X2, fit$labels)
  if (is.finite(ceiling)) {
    expected <- N$bands_profile(final_stats, ceiling, 1.1, 2)$J
  } else {
    stopifnot(all(final_stats$rank == 2))
    expected <- log(2 * pi) + 1 +
      sum(final_stats$mass * final_stats$logdet) / 2 -
      1.1 * sum(final_stats$mass * log(final_stats$mass))
  }
  near(tail(fit$trace, 1), expected, 1e-9)
}

# The profiled finite-C criterion must be nonincreasing as C increases.
scores <- vapply(c(.01, .1, 1, 10, 100), function(C) {
  N$bands_profile(a, C, 1.1, 2)$J
}, 0)
stopifnot(all(diff(scores) <= 1e-10))
