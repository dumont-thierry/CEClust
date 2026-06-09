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

test_that("CECclassif skips non-finite shots in mixed base mode", {
  Z <- simulate_multidim_benchmark_data(n = 500, p_num = 4, p_fac = 2, seed = 11)

  CECconfigure_runtime("base")
  set.seed(20260525L)
  fit <- NULL
  expect_no_error(
    fit <- CECclassif(
      Z = Z,
      lambda = 1,
      C = 10,
      r0 = 10,
      Nshots = 2,
      Nloop = 5,
      familyType = "gaussAndDiscreteVector",
      backend_data = list(optimized = FALSE, raw = Z)
    )
  )

  expect_true(is.null(fit) || is.finite(fit$Hphi))
  if (!is.null(fit)) {
    expect_equal(length(fit$clusters), nrow(Z))
  }
})

test_that("independent shot strategy is reproducible across core counts for non-1D families", {
  CECconfigure_runtime("fast")

  set.seed(11)
  cases <- list(
    gaussVector = list(
      Z = rbind(
        matrix(rnorm(40, -2), ncol = 2),
        matrix(rnorm(40, 2), ncol = 2)
      ),
      familyType = "gaussVector"
    ),
    discreteVector = list(
      Z = data.frame(
        a = factor(rep(c("a", "b"), each = 20)),
        b = factor(rep(c("x", "y"), each = 20))
      ),
      familyType = "discreteVector"
    ),
    gaussAndDiscreteVector = list(
      Z = data.frame(
        x = c(rnorm(20, -2), rnorm(20, 2)),
        y = c(rnorm(20, 0), rnorm(20, 3)),
        f = factor(rep(c("a", "b"), each = 20))
      ),
      familyType = "gaussAndDiscreteVector"
    )
  )

  for (case in cases) {
    set.seed(321)
    sequential <- CECclassif(
      Z = case$Z,
      familyType = case$familyType,
      C = 10,
      r0 = 4,
      Nshots = 3,
      Nloop = 8,
      shot_strategy = "independent",
      n_cores = 1
    )

    set.seed(321)
    parallel <- CECclassif(
      Z = case$Z,
      familyType = case$familyType,
      C = 10,
      r0 = 4,
      Nshots = 3,
      Nloop = 8,
      shot_strategy = "independent",
      n_cores = 2
    )

    expect_true(is.finite(sequential$Hphi))
    expect_equal(parallel$Hphi, sequential$Hphi, tolerance = 1e-12)
    expect_equal(parallel$REO, sequential$REO)
    expect_equal(parallel$clusters, sequential$clusters)
  }
})

test_that("CECclassif returns coherent public fits for non-1D families", {
  CECconfigure_runtime("fast")

  set.seed(12)
  cases <- list(
    gaussVector = list(
      Z = rbind(
        matrix(rnorm(40, -2), ncol = 2),
        matrix(rnorm(40, 2), ncol = 2)
      ),
      familyType = "gaussVector"
    ),
    discreteVector = list(
      Z = data.frame(
        a = factor(rep(c("a", "b"), each = 20)),
        b = factor(rep(c("x", "y"), each = 20))
      ),
      familyType = "discreteVector"
    ),
    gaussAndDiscreteVector = list(
      Z = data.frame(
        x = c(rnorm(20, -2), rnorm(20, 2)),
        y = c(rnorm(20, 0), rnorm(20, 3)),
        f = factor(rep(c("a", "b"), each = 20))
      ),
      familyType = "gaussAndDiscreteVector"
    )
  )

  for (case_name in names(cases)) {
    case <- cases[[case_name]]
    n <- NROW(case$Z)

    set.seed(412)
    fit <- CECclassif(
      Z = case$Z,
      familyType = case$familyType,
      C = 10,
      r0 = 4,
      Nshots = 2,
      Nloop = 8
    )

    expect_true(is.finite(fit$Hphi), info = case_name)
    expect_equal(length(fit$clusters), n, info = case_name)
    expect_equal(length(fit$phi), n * fit$REO, info = case_name)
    expect_equal(length(fit$params$phi), length(fit$phi), info = case_name)
    expect_equal(unclass(fit$params$phi), unclass(fit$phi), info = case_name)
    expect_equal(fit$params$familyType, case$familyType, info = case_name)
    expect_equal(length(fit$params$states), fit$REO, info = case_name)
    expect_equal(length(fit$params$nu), fit$REO, info = case_name)
    expect_equal(sum(fit$params$nu), 1, tolerance = 1e-12, info = case_name)
    expect_true(all(fit$clusters >= 1L & fit$clusters <= fit$REO), info = case_name)

    if (case$familyType == "gaussVector") {
      expect_equal(nrow(fit$params$m), fit$REO)
      expect_equal(ncol(fit$params$m), NCOL(case$Z))
      expect_equal(length(fit$params$Sigma), fit$REO)
    }

    if (case$familyType == "discreteVector") {
      expect_identical(fit$params$Cquali, Inf)
      expect_equal(length(fit$params$discreteProbList), fit$REO)
      expect_equal(fit$params$qualitativeBoundReached, integer(0))
    }

    if (case$familyType == "gaussAndDiscreteVector") {
      expect_identical(fit$params$Cquali, Inf)
      expect_equal(unname(fit$params$colNum), 1:2)
      expect_equal(unname(fit$params$colFactor), 3L)
      expect_equal(length(fit$params$Sigma), fit$REO)
      expect_equal(length(fit$params$discreteProbList), fit$REO)
      expect_equal(fit$params$qualitativeBoundReached, integer(0))
    }
  }
})

test_that("one-shot entropy plotting remains compatible with legacy and fast paths", {
  CECconfigure_runtime("fast")

  tmp <- tempfile(fileext = ".pdf")
  grDevices::pdf(tmp)
  on.exit({
    if (grDevices::dev.cur() > 1L) {
      grDevices::dev.off()
    }
    unlink(tmp)
  }, add = TRUE)

  set.seed(31)
  Z_univ <- c(rnorm(10, -1), rnorm(10, 1))
  fit_univ <- NULL
  expect_no_error(
    fit_univ <- CECclassifOneShot(
      Z = Z_univ,
      familyType = "gaussUniv",
      r0 = 3,
      Nloop = 4,
      displayPlotEntropy = TRUE
    )
  )
  expect_true(is.finite(fit_univ$Hphi))
  expect_equal(length(fit_univ$phi) %% length(Z_univ), 0)

  set.seed(32)
  Z_mixed <- data.frame(
    x = c(rnorm(10, -1), rnorm(10, 1)),
    f = factor(rep(c("a", "b"), each = 10))
  )
  fit_mixed <- NULL
  expect_no_error(
    fit_mixed <- CECclassifOneShot(
      Z = Z_mixed,
      familyType = "gaussAndDiscreteVector",
      C = 10,
      r0 = 3,
      Nloop = 4,
      displayPlotEntropy = TRUE
    )
  )
  expect_true(is.finite(fit_mixed$Hphi))
  expect_equal(length(fit_mixed$phi) %% nrow(Z_mixed), 0)
})

test_that("linked-lambda workflow returns extractable best partitions", {
  Z <- simulate_multidim_benchmark_data(n = 60, p_num = 3, p_fac = 2, seed = 11)

  lambda_diag <- CECfitLambdaGrid(
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

test_that("stable-lambda selection can exclude saturation-heavy lambdas", {
  lambda_diag <- list(
    summary = data.frame(
      lambda = c(0.1, 0.2, 0.3),
      stab_ratio_mean = c(0.70, 0.95, 0.95),
      stabB_mean = c(0.9, 0.9, 0.9)
    )
  )

  best_parts <- list(
    summary = data.frame(lambda = c(0.1, 0.2, 0.3)),
    best = list(
      list(params = list(nu = c(0.85, 0.15), functionBoundReached = 1L)),
      list(params = list(nu = c(0.20, 0.80), functionBoundReached = 1L)),
      list(params = list(nu = c(0.10, 0.90), functionBoundReached = integer(0)))
    )
  )

  stable <- CECidentifyStableLambdas(
    lambda_diag,
    rule = "both",
    ratio_threshold = 0.8,
    stabB_threshold = 0.8,
    sat_threshold = 0.8,
    use_smoothed = FALSE,
    best_partitions_obj = best_parts
  )

  expect_equal(stable$summary$sat, c(0.85, 0.2, 0))
  expect_equal(stable$summary$keep_stability, c(FALSE, TRUE, TRUE))
  expect_equal(stable$summary$exclude_saturation, c(TRUE, FALSE, FALSE))
  expect_equal(stable$summary$stable, c(FALSE, TRUE, TRUE))
  expect_equal(stable$lambda_min, 0.2)
})

test_that("PCA trajectory helpers work for multivariate numeric data", {
  skip_if_not_installed("ggplot2")

  Z <- iris[, -5]

  lambda_diag <- CECfitLambdaGrid(
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

test_that("CECpredict predicts all missing numeric coordinates cluster-wise", {
  CECpredict_internal <- getFromNamespace("CECpredict", "CEClust")

  params <- list(
    familyType = "gaussAndDiscreteVector",
    states = 1:2,
    nu = c(0.5, 0.5),
    lambda = 1,
    C = 10,
    colFactor = 1,
    colNum = 2:3,
    factors = c("a", "b"),
    discreteProbList = list(
      matrix(c(0.9, 0.1), nrow = 2, ncol = 1),
      matrix(c(0.1, 0.9), nrow = 2, ncol = 1)
    ),
    m = matrix(c(10, 100, 20, 200), nrow = 2, byrow = TRUE),
    Sigma = list(diag(2), diag(2))
  )

  Zpred <- data.frame(f = factor(c("a", "b", "a"), levels = c("a", "b")))
  pred <- CECpredict_internal(Zpred, params, idColToPred = c(2, 3))

  expect_equal(unname(as.matrix(pred)), matrix(c(10, 100, 20, 200, 10, 100), ncol = 2, byrow = TRUE))
})
