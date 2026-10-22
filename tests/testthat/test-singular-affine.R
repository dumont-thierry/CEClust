test_that("affine scores include new points in the translated span", {
  fit <- getFromNamespace("adjust_covariance", "CEClust")
  density <- getFromNamespace("CECgaussian_reference_logdens", "CEClust")
  x <- seq(-2, 2, length.out = 9)
  X <- cbind(x, 2*x + 3)
  queries <- rbind(X, c(100, 203), c(0, 3.01))
  for (C in c(.01, .5, 1, 10)) {
    b <- rep(1/(sqrt(2*pi)*sqrt(C)), 2)
    S <- fit(X, b)$Sigma_corrected
    values <- exp(density(queries, colMeans(X), S))
    expect_equal(values[1:10], rep(C, 10), tolerance=1e-11)
    expect_equal(values[11], 0, tolerance=1e-11)
    expect_equal(as.numeric(S), as.numeric(cov(X) * (nrow(X)-1)/nrow(X)))
  }
})

test_that("zero covariance uses the singleton affine support", {
  f <- getFromNamespace("adjust_covariance", "CEClust")
  density <- getFromNamespace("CECgaussian_reference_logdens", "CEClust")
  X <- matrix(rep(c(5, -7, 2), each=6), ncol=3)
  for (C in c(.2, 4)) {
    S <- f(X, rep(1/(sqrt(2*pi)*C^(1/3)), 3))$Sigma_corrected
    values <- exp(density(rbind(X[1,], c(5,-7,2.1)), X[1,], S))
    expect_equal(values, c(C, 0))
    expect_equal(attr(S, "CEC_singular_support")$rank, 0L)
  }
})

test_that("base and compiled scores agree for rotated affine supports", {
  set.seed(8801)
  f <- getFromNamespace("CECoptParam_fast_from_clusters", "CEClust")
  density <- getFromNamespace("CECcompute_logdens_matrix_fast", "CEClust")
  toggle <- getFromNamespace("CECset_fast_backend", "CEClust")
  on.exit(toggle(TRUE))
  for (d in c(2, 8)) {
    Q <- qr.Q(qr(matrix(rnorm(d*d), d)))
    X <- outer(seq(-3,3,length.out=15), Q[,1])
    offset <- seq_len(d)*100
    X <- sweep(X, 2, offset, "+")
    queries <- rbind(X, offset + 10*Q[,1], offset + .05*Q[,d])
    backend <- CECprepare_backend_data(X, "gaussVector")
    for (C in c(.3, 7)) {
      toggle(FALSE)
      params <- f(backend, rep(1L,nrow(X)), C=C)
      query_backend <- CECprepare_backend_data(queries, "gaussVector")
      base <- density(query_backend, params)
      toggle(TRUE)
      fast <- density(query_backend, params)
      expect_equal(fast, base, tolerance=1e-11)
      expect_equal(as.numeric(exp(base[1:16,])), rep(C,16), tolerance=1e-10)
      expect_equal(exp(base[17,1]), 0)
    }
  }
})

test_that("entropy ignores zero weights at zero off-support scores", {
  x <- seq(-2,2,length.out=7)
  X <- rbind(cbind(x, x*0), cbind(x, x*0+4))
  labels <- rep(1:2, each=7)
  phi <- getFromNamespace("CECclusters_to_phi", "CEClust")(labels, 2)
  evaluate <- getFromNamespace("evalCompositeEntropy", "CEClust")
  for (lambda in c(.7, 1, 2)) {
    expect_equal(evaluate(phi, X, lambda=lambda, C=.2, familyType="gaussVector"),
                 log(2)-log(.2)/lambda, tolerance=1e-10)
  }
})

test_that("mixed-data and bootstrap projection paths use the same scores", {
  x <- as.numeric(seq(-2,2,length.out=7))
  X <- data.frame(x=rep(x,2), y=rep(c(0,4), each=7),
                  q=factor(rep(c("a","b"), each=7)))
  labels <- rep(1:2,each=7)
  phi <- getFromNamespace("CECclusters_to_phi", "CEClust")(labels, 2)
  backend <- CECprepare_backend_data(X, "gaussAndDiscreteVector")
  f <- getFromNamespace("CECoptParam_fast_from_clusters", "CEClust")
  params <- f(backend, labels, C=.2)
  fast <- getFromNamespace("CECcompute_logdens_matrix_fast", "CEClust")(backend, params)
  base <- matrix(getFromNamespace("CECdens", "CEClust")(X, params, applyLog=TRUE)$density, nrow(X))
  expect_equal(fast, base)
  expect_equal(as.numeric(exp(fast[cbind(seq_len(nrow(X)),labels)])), rep(.2, nrow(X)))
  expect_equal(params$functionBoundReached, 1:2)
  projected <- getFromNamespace("CECproject_fit_on_reference", "CEClust")(
    list(params=params, phi=phi), X)
  expect_equal(projected$clusters, labels)
  expect_equal(projected$H, log(2)-log(.2), tolerance=1e-10)
})

test_that("new observations outside singular supports use the proposal floor", {
  X <- cbind(as.numeric(-2:2), rep(3,5))
  backend <- CECprepare_backend_data(X, "gaussVector")
  params <- getFromNamespace("CECoptParam_fast_from_clusters", "CEClust")(
    backend, rep(1L,5), C=.2)
  query <- matrix(c(0,4),nrow=1)
  expect_equal(as.numeric(getFromNamespace("optPhi", "CEClust")(query, params)), 1)
  step <- getFromNamespace("CECoptPhi_backend", "CEClust")(
    query, params, backend_data=CECprepare_backend_data(query,"gaussVector"))
  expect_equal(step$clusters, 1L)
  expect_equal(as.numeric(step$assigned_logdens), -Inf)
  projected <- getFromNamespace("CECproject_fit_on_reference", "CEClust")(
    list(params=params), query)
  expect_equal(projected$clusters, 1L)
  expect_equal(projected$H, -log(.2))
})

test_that("singular metadata survives serialization and public fitting", {
  X <- cbind(as.numeric(-4:4), rep(3,9))
  for (C in c(.2, 10)) for (runtime in c("base", "fast")) {
    CECconfigure_runtime(runtime)
    set.seed(812)
    result <- CECclassif(X, C=C, r0=2, Nshots=2, Nloop=10, familyType="gaussVector")
    expect_true(is.finite(result$Hphi))
    expect_equal(result$determinant_diagnostics$rejected_shots, 0L)
    expect_true(length(result$singular_clusters) > 0)
    expect_equal(result$params$functionBoundReached, seq_along(result$params$nu))
    expect_true(all(diff(result$descent_diagnostics$H_trace) <= 0))
    restored <- unserialize(serialize(result,NULL))
    expect_equal(getFromNamespace("optPhi", "CEClust")(X,restored$params),
                 getFromNamespace("optPhi", "CEClust")(X,result$params))
    H <- getFromNamespace("evalCompositeEntropy", "CEClust")(
      result$phi, X, lambda=1, C=C, familyType="gaussVector")
    expect_equal(result$Hphi,H,tolerance=1e-10)
  }
})

test_that("bootstrap projection remains defined outside fitted affine supports", {
  X <- rbind(c(0,0),c(1,0),c(0,1))
  grid <- CECfitLambdaGrid(X, lambda_grid=1, C=10, r0=1, k0=1, B=8,
    Nshots_fresh=1, Nshots_warm=1, Nloop=3, familyType="gaussVector",
    seed=451, silent=TRUE, verbose=FALSE, n_cores=1L,
    checkpoint_dir=FALSE, auto_checkpoint=FALSE, show_progress=FALSE)
  expect_equal(grid$summary$n_validB, 8)
  # A single original run cannot estimate exhaustive optimization stability.
  expect_true(is.na(grid$summary$stab_ratio_mean))
  expect_equal(grid$summary$stabB_mean, 1)
  fits <- grid$details[[1]]$fitsB
  singular <- which(vapply(fits, function(fit) !is.null(
    attr(fit$params$Sigma[[1]], "CEC_singular_support")), logical(1)))
  expect_true(length(singular) > 0)
  for (b in singular) {
    expect_equal(grid$details[[1]]$partitionsB[[b]], rep(1L,3))
    expect_true(is.finite(grid$details[[1]]$HB[b]))
  }
})
