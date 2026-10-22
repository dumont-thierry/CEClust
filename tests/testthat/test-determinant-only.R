fit_cov <- getFromNamespace("CECfit_covariance_determinant", "CEClust")
fit_weighted <- getFromNamespace("adjust_covariance", "CEClust")
fit_matrix <- getFromNamespace("CECadjust_covariance_from_matrix", "CEClust")
logdet <- function(S) 2 * sum(log(diag(chol(S))))
cost <- function(S, A) logdet(A) / 2 + sum(diag(solve(A, S))) / 2

test_that("one-dimensional fit, including small groups, is the clipped MLE", {
  for (n in c(2, 3, 10)) for (b in c(0.01, 0.7, 3)) {
    x <- matrix(seq_len(n), ncol = 1)
    v <- mean((x - mean(x))^2)
    got <- fit_weighted(x, b)
    expect_equal(as.numeric(got$Sigma_corrected), max(v, b^2), tolerance = 1e-12)
  }
  x <- matrix(c(0, 2, 6), ncol = 1)
  expect_equal(fit_weighted(x, 1, c(1, 2, 1))$Sigma_corrected,
               fit_weighted(x, 1, c(.25, .5, .25))$Sigma_corrected)
})

test_that("full-rank determinant projection is exact and rotation equivariant", {
  set.seed(147)
  for (d in c(2, 5, 10, 20)) {
    A <- matrix(rnorm(d * d), d)
    S <- crossprod(A) + diag(d)
    Q <- qr.Q(qr(matrix(rnorm(d * d), d)))
    for (relative_volume in c(.1, 1, 10, 10000)) {
      ldmin <- logdet(S) + log(relative_volume)
      b <- rep(exp(ldmin / (2 * d)), d)
      got <- fit_cov(S, b)
      expected <- S * max(1, relative_volume^(1/d))
      expect_equal(got$Sigma_corrected, expected, tolerance = 1e-10)
      expect_gte(logdet(got$Sigma_corrected) - ldmin, -1e-10)
      expect_equal(got$boundReached, relative_volume >= 1)
      rotated <- fit_cov(Q %*% S %*% t(Q), b)
      expect_equal(rotated$Sigma_corrected, Q %*% expected %*% t(Q), tolerance = 1e-9)
    }
  }
})

test_that("two-dimensional optimizer agrees with independent numerical minimization", {
  S <- matrix(c(.1, .04, .04, .08), 2)
  delta <- .25
  sigma <- function(z) {
    a <- exp(z[1])
    b <- z[2]
    matrix(c(a, b, b, (delta + b^2) / a), 2)
  }
  optimum <- optim(c(0, -.4), function(z) cost(S, sigma(z)),
                   method = "BFGS", control = list(reltol = 1e-12, maxit = 1000))
  exact <- fit_cov(S, rep(delta^.25, 2))$Sigma_corrected
  expect_equal(sigma(optimum$par), exact, tolerance = 1e-5)
  expect_equal(optimum$value, cost(S, exact), tolerance = 1e-10)
})

test_that("log determinants handle overflow and underflow of raw products", {
  for (exponent in c(-100, 100)) {
    d <- 40
    S <- diag(exp(2*exponent - 20), d)
    out <- fit_cov(S, rep(exp(exponent), d))
    expect_true(all(is.finite(out$Sigma_corrected)))
    expect_equal(out$log_det_fitted, 2 * exponent * d, tolerance = 1e-10)
  }
})

test_that("small full-rank groups use their covariance; no diagonal fallback", {
  X <- rbind(c(0, 0), c(2, 0), c(0, 1))
  S <- crossprod(sweep(X, 2, colMeans(X))) / nrow(X)
  out <- fit_weighted(X, c(.1, .1))
  expect_equal(out$Sigma_corrected, S)
  expect_equal(fit_matrix(S, c(.1, .1), count = 3)$Sigma_corrected, S)
})

test_that("zero covariance attains the cap even when C is small", {
  d <- 2
  C <- .01
  b <- rep(1/(sqrt(2*pi)*C^(1/d)), d)
  out <- fit_weighted(matrix(.1, 7, d), b)
  density <- getFromNamespace("CECgaussian_reference_logdens", "CEClust")
  peak <- exp(density(matrix(.1, 1, d), rep(.1, d), out$Sigma_corrected))
  expect_equal(as.numeric(out$Sigma_corrected), rep(0, d*d))
  expect_equal(peak, C, tolerance = 1e-12)
  expect_true(out$boundReached)
})

test_that("singular nonzero covariance is retained with support metadata", {
  for (S in list(diag(c(1,0)), diag(c(1,1e-20)))) {
    out <- fit_cov(S, c(1,1))
    expect_equal(as.numeric(out$Sigma_corrected), as.numeric(S))
    expect_equal(out$empirical_rank, 1L)
    expect_true(out$boundReached)
    expect_equal(out$covariance_status, "singular_affine_support")
  }
  expect_error(fit_cov(matrix(c(1,2,2,1), 2), c(1,1)), "positive semidefinite")
})

test_that("old rectangular calibration counterexample is corrected", {
  for (d in c(2, 10, 20)) {
    C <- 20 * (6*.99/pi)^(d/2)
    b <- rep(1/(sqrt(2*pi)*C^(1/d)), d)
    S <- diag(c(1/(12*20^2), rep(1/12, d-1)))
    got <- fit_matrix(S, b, count = 10000)
    peak <- exp(-d/2*log(2*pi) - logdet(got$Sigma_corrected)/2)
    expect_equal(peak, C, tolerance = 1e-10)
    expect_true(got$boundReached)
    expect_equal(got$Sigma_corrected, S / .99, tolerance = 1e-10)
  }
})

test_that("base and compiled fitting agree for the same fixed partition", {
  set.seed(161)
  X <- rbind(matrix(rnorm(400), ncol = 2),
             matrix(rnorm(400), ncol = 2) + 4)
  labels <- rep(1:2, each = 200)
  prepared <- getFromNamespace("CECprepare_backend_data", "CEClust")(X, "gaussVector")
  toggle <- getFromNamespace("CECset_fast_backend", "CEClust")
  fit_params <- getFromNamespace("CECoptParam_fast_from_clusters", "CEClust")
  toggle(FALSE)
  base <- fit_params(prepared, labels, C = .02)
  toggle(TRUE)
  fast <- fit_params(prepared, labels, C = .02)
  expect_equal(base$Sigma, fast$Sigma, tolerance = 1e-11)
  expect_equal(base$functionBoundReached, fast$functionBoundReached)
  for (A in fast$Sigma) {
    expect_equal(exp(-log(2*pi)-logdet(A)/2), .02, tolerance = 1e-10)
  }
})

test_that("public CECclassif produces finite saturated and unsaturated fits", {
  set.seed(200)
  X <- matrix(rnorm(1200), ncol = 2)
  for (C in c(.01, 1)) {
    out <- CECclassif(X, lambda = 1, C = C, r0 = 1, Nshots = 1,
                     Nloop = 10, familyType = "gaussVector")
    expect_true(is.finite(out$Hphi))
    expect_equal(out$REO, 1L)
    expect_equal(length(out$params$functionBoundReached), as.integer(C == .01))
  }
})

test_that("compiled and base code preserve zero covariance for identical points", {
  X <- matrix(.1, nrow = 7, ncol = 2)
  labels <- rep(1L, nrow(X))
  prepared <- getFromNamespace("CECprepare_backend_data", "CEClust")(X, "gaussVector")
  toggle <- getFromNamespace("CECset_fast_backend", "CEClust")
  f <- getFromNamespace("CECoptParam_fast_from_clusters", "CEClust")
  for (enabled in c(FALSE, TRUE)) {
    toggle(enabled)
    params <- f(prepared, labels, C = .01)
    expect_equal(params$functionBoundReached, 1L)
    density <- getFromNamespace("CECcompute_logdens_matrix_fast", "CEClust")
    expect_equal(as.numeric(exp(density(prepared, params))), rep(.01, nrow(X)), tolerance=1e-11)
    expect_equal(as.numeric(params$Sigma[[1]]), rep(0,4))
  }
})

test_that("singular shots are retained as saturated fits", {
  X <- rbind(c(-1, 0), c(0, 0), c(1, 0))
  shot <- getFromNamespace("CECrun_one_shot_with_regroup", "CEClust")(
    X, C=1, r0=1, Nloop=2, familyType="gaussVector")
  expect_equal(shot$Hphi, 0, tolerance=1e-12)
  expect_null(shot$determinant_fit_issue)
  out <- CECclassif(X, C=1, r0=1, Nshots=1, Nloop=2, familyType="gaussVector")
  expect_equal(out$REO, 1L)
  expect_equal(out$singular_clusters, 1L)
  expect_equal(out$determinant_diagnostics$rejected_shots, 0L)
  set.seed(11)
  X <- rbind(matrix(rnorm(40, -2), ncol=2), matrix(rnorm(40, 2), ncol=2))
  set.seed(412)
  fit <- CECclassif(X, C=10, r0=4, Nshots=3, Nloop=8, familyType="gaussVector")
  expect_true(is.finite(fit$Hphi))
  expect_equal(fit$determinant_diagnostics$attempted_shots, 3L)
  expect_gte(fit$determinant_diagnostics$rejected_shots, 0L)
})

test_that("lambda one calibration gives the largest rectangular minimizer", {
  for (d in c(1,2,10,20)) {
    C0 <- (6/pi)^(d/2)
    C <- 20*C0
    F <- function(r) ifelse(r<=20, -log(C0)+d/2,
                            log(r)-log(C)+d/2*(20/r)^(2/d))
    expect_equal(F(1:20), rep(F(1),20))
    expect_gt(F(21), F(20))
  }
})
