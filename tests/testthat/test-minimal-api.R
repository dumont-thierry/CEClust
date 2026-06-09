test_that("minimal API entry points are available", {
  expect_true(is.function(CECfitBoundGrid))
  expect_true(is.function(CECfitPreset))
  expect_true(is.function(CECselectStableLambdas))
  expect_true(is.function(CECsummariseGrid))
  expect_true(is.function(CECplotGrid))
  expect_true(is.function(CECplotPartition))
  expect_true(is.function(CECplotPath))
  expect_true(is.function(CECexplore))
})

test_that("preset names are intentionally small", {
  expect_equal(
    eval(formals(CECfitPreset)$preset),
    c("gaussian_1d", "iris", "breast_cancer")
  )
})

test_that("CECplotPath exposes diagnostic and partition views", {
  expect_equal(eval(formals(CECplotPath)$type), c("summary", "partitions"))
  expect_true(all(c(
    "stab_algo_threshold",
    "rsi_threshold",
    "sat_threshold",
    "lambda_stable_obj"
  ) %in% names(formals(CECplotPath))))
})

test_that("public grid wrappers do not write files by default", {
  bound_defaults <- as.list(formals(CECfitBoundGrid))
  preset_defaults <- as.list(formals(CECfitPreset))

  for (defaults in list(bound_defaults, preset_defaults)) {
    expect_false(defaults$save_results)
    expect_false(defaults$checkpoint_dir)
    expect_false(defaults$auto_checkpoint)
    expect_false(defaults$resume)
    expect_false(defaults$force_recompute)
  }
})

test_that("qualitative floors are separate from density saturation", {
  expect_identical(eval(formals(CECclassif)$Cquali), Inf)
  expect_identical(eval(formals(CECfitLambdaGrid)$Cquali), Inf)

  Z <- data.frame(q = factor(c("a", "b", "b", "b", "b")))
  backend <- CECprepare_backend_data(Z, familyType = "discreteVector")
  clusters <- rep(1L, nrow(Z))

  params_inf <- CECoptParam_fast_from_clusters(
    backend,
    clusters,
    C = 1,
    Cquali = Inf
  )
  expect_equal(params_inf$functionBoundReached, integer(0))
  expect_equal(params_inf$qualitativeBoundReached, integer(0))
  expect_equal(sort(params_inf$discreteProbList[[1]][, 1]), c(0.2, 0.8))

  params_floor <- CECoptParam_fast_from_clusters(
    backend,
    clusters,
    C = 1,
    Cquali = 3
  )
  expect_equal(params_floor$functionBoundReached, integer(0))
  expect_equal(params_floor$qualitativeBoundReached, 1L)
  expect_equal(
    sort(params_floor$discreteProbList[[1]][, 1]),
    c(1 / 3, 2 / 3),
    tolerance = 1e-12
  )
})

test_that("mixed fast params keep qualitative and gaussian bounds separate", {
  Z_qual <- data.frame(
    x = c(0, 0, 0, 0, 0),
    q = factor(c("a", "b", "b", "b", "b"))
  )
  backend_qual <- CECprepare_backend_data(Z_qual, familyType = "gaussAndDiscreteVector")
  clusters_qual <- rep(1L, nrow(Z_qual))

  params_inf <- CECoptParam_fast_from_clusters(
    backend_qual,
    clusters_qual,
    C = 10,
    Cquali = Inf
  )
  expect_identical(params_inf$Cquali, Inf)
  expect_equal(params_inf$functionBoundReached, integer(0))
  expect_equal(params_inf$qualitativeBoundReached, integer(0))
  expect_equal(sort(params_inf$discreteProbList[[1]][, 1]), c(0.2, 0.8))

  params_floor <- CECoptParam_fast_from_clusters(
    backend_qual,
    clusters_qual,
    C = 10,
    Cquali = 3
  )
  expect_identical(params_floor$Cquali, 3)
  expect_equal(params_floor$functionBoundReached, integer(0))
  expect_equal(params_floor$qualitativeBoundReached, 1L)
  expect_equal(
    sort(params_floor$discreteProbList[[1]][, 1]),
    c(1 / 3, 2 / 3),
    tolerance = 1e-12
  )

  Z_small <- data.frame(
    x = c(0, 0, 1, 1, 1),
    q = factor(c("a", "b", "b", "b", "b"))
  )
  backend_small <- CECprepare_backend_data(Z_small, familyType = "gaussAndDiscreteVector")
  params_small <- CECoptParam_fast_from_clusters(
    backend_small,
    clusters = c(1L, 1L, 2L, 2L, 2L),
    C = 10,
    Cquali = Inf
  )

  expect_identical(params_small$Cquali, Inf)
  expect_equal(params_small$functionBoundReached, 1:2)
  expect_equal(params_small$qualitativeBoundReached, integer(0))
  expect_equal(unname(params_small$colNum), 1L)
  expect_equal(unname(params_small$colFactor), 2L)
})

test_that("C++ cluster choice from log densities preserves assigned log densities", {
  set.seed(91)
  logdens <- matrix(rnorm(30), nrow = 10, ncol = 3)
  nu <- c(0.2, 0.3, 0.5)
  lambda <- 0.7

  score_choice <- cec_cpp_choose_clusters(logdens / lambda, nu)
  logdens_choice <- cec_cpp_choose_clusters_from_logdens(logdens, nu, lambda)

  expect_equal(logdens_choice$clusters, score_choice$clusters)
  expect_equal(
    logdens_choice$assigned_logdens,
    logdens[cbind(seq_len(nrow(logdens)), logdens_choice$clusters)],
    tolerance = 1e-12
  )
})
