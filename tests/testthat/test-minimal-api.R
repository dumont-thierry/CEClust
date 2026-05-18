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
