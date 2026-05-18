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
