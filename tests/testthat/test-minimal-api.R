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

