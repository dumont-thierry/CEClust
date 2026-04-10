test_that("runtime configuration returns the expected structure", {
  runtime <- CECconfigure_runtime("base")

  expect_s3_class(runtime, "CEC_runtime_config")
  expect_named(
    runtime,
    c("profile", "use_cpp", "n_cores", "max_cores", "prefer_multicore", "performance")
  )
  expect_equal(runtime$profile, "base")
  expect_equal(runtime$n_cores, 1L)
})

test_that("linked-lambda workflow returns extractable best partitions", {
  Z <- simulate_multidim_benchmark_data(n = 60, p_num = 3, p_fac = 2, seed = 11)

  lambda_diag <- CECdiagnose_lambda_grid_linked(
    Z = Z,
    lambda_grid = seq(0.2, 0.6, by = 0.2),
    k0 = 2,
    B = 2,
    C = 10,
    r0 = 6,
    Nshots_fresh = 1,
    Nshots_warm = 1,
    Nloop = 10,
    familyType = "gaussAndDiscreteVector",
    seed = 11,
    silent = TRUE,
    verbose = FALSE,
    n_cores = 1L,
    checkpoint_dir = FALSE,
    auto_checkpoint = FALSE,
    show_progress = FALSE
  )

  expect_true(all(c("summary", "details", "lambda_grid", "linked", "execution") %in% names(lambda_diag)))

  lambda_diag <- CECaddSmoothedDiagnostics(lambda_diag, k = 3)
  stable <- CECidentifyStableLambdas(lambda_diag)
  best_parts <- CECextractBestPartitions(lambda_diag, source = "fixed", Z = Z)
  coherence <- CECcheckBestPartitionCoherence(best_parts, Z)
  changes <- CECdetectPartitionChanges(best_parts)
  repaired <- CECrepairBestPartitionTrajectory(best_parts, Z, max_iter = 10)

  expect_true("summary" %in% names(stable))
  expect_true(all(c("best", "summary", "source", "criterion") %in% names(best_parts)))
  expect_true(is.data.frame(best_parts$summary))
  expect_true("all_coherent" %in% names(coherence))
  expect_true("change_lambdas" %in% names(changes))
  expect_true(all(c("repair_log", "coherence_before", "coherence_after") %in% names(repaired)))
})

test_that("PCA trajectory helpers work for multivariate numeric data", {
  Z <- iris[, -5]

  lambda_diag <- CECdiagnose_lambda_grid_linked(
    Z = Z,
    lambda_grid = seq(0.2, 0.4, by = 0.2),
    k0 = 2,
    B = 1,
    C = 10,
    r0 = 4,
    Nshots_fresh = 1,
    Nshots_warm = 1,
    Nloop = 10,
    familyType = "gaussVector",
    seed = 5,
    silent = TRUE,
    verbose = FALSE,
    n_cores = 1L,
    checkpoint_dir = FALSE,
    auto_checkpoint = FALSE,
    show_progress = FALSE
  )

  best_parts <- CECextractBestPartitions(lambda_diag, source = "fixed", Z = Z)
  path_df <- build_partition_path_df(Z, best_parts)
  p <- plot_partition_path_pca(path_df)

  expect_true(all(c("PC1", "PC2", "id", "lambda", "cluster", "REO") %in% names(path_df)))
  expect_s3_class(p, "ggplot")
})
