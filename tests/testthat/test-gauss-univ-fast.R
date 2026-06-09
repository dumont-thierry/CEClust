test_that("gaussUniv returns a valid public fit with fast runtime enabled", {
  CECconfigure_runtime("fast")

  set.seed(101)
  Z <- c(rnorm(30, -2), rnorm(30, 2))

  fit <- CECclassif(
    Z = Z,
    familyType = "gaussUniv",
    r0 = 5,
    Nshots = 2,
    Nloop = 10
  )

  expect_true(is.finite(fit$Hphi))
  expect_equal(length(fit$clusters), length(Z))
  expect_equal(length(fit$phi), length(Z) * fit$REO)
  expect_equal(length(fit$params$states), fit$REO)
  expect_equal(length(fit$params$m), fit$REO)
  expect_equal(length(fit$params$s), fit$REO)
  expect_equal(length(fit$params$phi), length(fit$phi))
  expect_equal(unclass(fit$params$phi), unclass(fit$phi))
  expect_equal(
    attr(fit$params$phi, "CEC_clusters", exact = TRUE),
    attr(fit$phi, "CEC_clusters", exact = TRUE)
  )
  expect_equal(
    fit$clusters,
    attr(fit$phi, "CEC_clusters", exact = TRUE)
  )
  expect_true(all(fit$params$s >= 1 / (sqrt(2 * pi) * fit$params$C)))
  expect_equal(fit$params$familyType, "gaussUniv")
})

test_that("gaussUniv keeps legacy-compatible quality with fast runtime enabled", {
  set.seed(102)
  Z <- c(rnorm(40, -2), rnorm(40, 2))

  run_fit <- function(backend_data = NULL) {
    set.seed(202)
    CECclassif(
      Z = Z,
      familyType = "gaussUniv",
      r0 = 6,
      Nshots = 4,
      Nloop = 30,
      backend_data = backend_data
    )
  }

  legacy <- run_fit(list(optimized = FALSE, raw = Z))
  CECconfigure_runtime("fast")
  fast <- run_fit(NULL)

  expect_true(is.finite(legacy$Hphi))
  expect_true(is.finite(fast$Hphi))
  expect_lte(fast$Hphi, legacy$Hphi + 0.05)
})

test_that("gaussUniv independent shots are reproducible across core counts", {
  CECconfigure_runtime("fast")

  set.seed(112)
  Z <- c(rnorm(50, -2), rnorm(50, 2))

  set.seed(212)
  sequential <- CECclassif(
    Z = Z,
    familyType = "gaussUniv",
    r0 = 6,
    Nshots = 4,
    Nloop = 20,
    shot_strategy = "independent",
    n_cores = 1
  )

  set.seed(212)
  parallel <- CECclassif(
    Z = Z,
    familyType = "gaussUniv",
    r0 = 6,
    Nshots = 4,
    Nloop = 20,
    shot_strategy = "independent",
    n_cores = 2
  )

  expect_equal(parallel$Hphi, sequential$Hphi, tolerance = 1e-12)
  expect_equal(parallel$REO, sequential$REO)
  expect_equal(parallel$clusters, sequential$clusters)
})

test_that("gaussUniv fast entropy matches decomposed reference calculation", {
  set.seed(103)
  Z <- c(rnorm(50, -2), rnorm(50, 2))

  fit <- CECclassif(
    Z = Z,
    familyType = "gaussUniv",
    r0 = 6,
    Nshots = 3,
    Nloop = 20
  )

  H_ref <- evalCompositeEntropy(
    phi = fit$phi,
    Z = Z,
    lambda = fit$params$lambda,
    C = fit$params$C,
    familyType = "gaussUniv"
  )

  expect_equal(fit$Hphi, H_ref, tolerance = 1e-12)
})

test_that("gaussUniv final params match recomputation from final clusters", {
  set.seed(108)
  Z <- c(rnorm(45, -2), rnorm(45, 2), rnorm(30, 5))
  Z_num <- as.numeric(Z)
  Z_stats <- cbind(Z_num, Z_num * Z_num)

  CECconfigure_runtime("fast")
  set.seed(208)
  fit <- CECclassif(
    Z = Z,
    familyType = "gaussUniv",
    r0 = 7,
    Nshots = 4,
    Nloop = 30
  )

  recomputed <- CECoptParam_gaussUniv_from_clusters(
    Z = Z,
    clusters = fit$clusters,
    lambda = fit$params$lambda,
    C = fit$params$C,
    build_phi = TRUE,
    compress = TRUE,
    Z_stats = Z_stats
  )

  expect_equal(fit$params$nu, recomputed$nu, tolerance = 1e-12)
  expect_equal(fit$params$m, recomputed$m, tolerance = 1e-12)
  expect_equal(fit$params$s, recomputed$s, tolerance = 1e-12)
  expect_equal(fit$params$varPhi, recomputed$varPhi, tolerance = 1e-12)
  expect_equal(
    CECentropy_gaussUniv_from_params(fit$params),
    CECentropy_gaussUniv_from_params(recomputed),
    tolerance = 1e-12
  )
  expect_equal(
    attr(fit$params$phi, "CEC_clusters", exact = TRUE),
    fit$clusters
  )
})

test_that("C++ gaussUniv assignment matches matrix reference", {
  CECconfigure_runtime("fast")

  Z <- c(-3, -1, 0, 1, 2, 4)
  param <- list(
    familyType = "gaussUniv",
    nu = c(0.25, 0.5, 0.25),
    m = c(-2, 0.5, 3),
    s = c(0.7, 1.1, 0.9)
  )
  lambda <- 0.8
  scoreM <- numeric(length(Z) * length(param$nu))
  dim(scoreM) <- c(length(Z), length(param$nu))
  score_const <- log(param$nu) - log(param$s) / lambda
  score_quad <- -0.5 / (lambda * param$s * param$s)
  for (i in seq_along(param$nu)) {
    scoreM[, i] <- score_const[i] + score_quad[i] * (Z - param$m[i])^2
  }
  expected <- max.col(scoreM, ties.method = "first")

  expect_equal(
    cec_cpp_gauss_univ_assign(Z, param$nu, param$m, param$s, lambda),
    expected
  )
})

test_that("C++ gaussUniv random phi expansion matches R construction", {
  CECconfigure_runtime("fast")

  n <- 4L
  phi <- seq(0.1, 1.2, length.out = 12L)

  set.seed(777)
  expected <- c(phi * runif(length(phi)), runif(n))
  dim(expected) <- c(n, length(phi) / n + 1L)

  set.seed(777)
  actual <- cec_cpp_gauss_univ_expand_phi_random(phi, n)

  expect_equal(actual, expected, tolerance = 1e-15)
})

test_that("C++ gaussUniv warm perturbation matches R construction", {
  CECconfigure_runtime("fast")

  n <- 5L
  r <- 4L
  phi <- seq(0.02, 0.98, length.out = n * r)
  nu <- c(0.05, 0.4, 0.45, 0.1)

  set.seed(778)
  phiM <- matrix(phi, n, r)
  toKeep <- which(nu * n >= 10)
  if (length(toKeep) == 0) {
    toKeep <- which.max(nu)
  }
  expected <- phiM[, toKeep, drop = FALSE]
  expected <- expected + 0.1 * runif(length(expected))
  sums <- rowSums(expected)
  sums[sums <= 0] <- 1
  expected <- expected / sums
  next_expected <- runif(3)

  set.seed(778)
  actual <- cec_cpp_gauss_univ_perturb_warm_phiM(phi, r, nu, n)
  next_actual <- runif(3)

  expect_equal(actual, expected, tolerance = 1e-15)
  expect_equal(next_actual, next_expected, tolerance = 0)
})

test_that("C++ replacement sampler matches sample.int with replacement", {
  CECconfigure_runtime("fast")

  for (r in c(1L, 2L, 3L, 5L, 24L)) {
    for (n in c(1L, 10L, 5000L)) {
      set.seed(779)
      expected <- sample.int(r, n, replace = TRUE)
      next_expected <- runif(3)

      set.seed(779)
      actual <- cec_cpp_sample_int_replace(r, n)
      next_actual <- runif(3)

      expect_equal(actual, expected, info = paste("r", r, "n", n))
      expect_equal(next_actual, next_expected, tolerance = 0, info = paste("rng", r, n))
    }
  }
})

test_that("C++ gaussUniv cluster sums match rowsum reference", {
  CECconfigure_runtime("fast")

  Z <- c(-2, -1, 0, 1, 3, 5)
  clusters <- c(1L, 2L, 1L, 3L, 2L, 3L)
  r <- 3L
  stats <- cec_cpp_gauss_univ_cluster_sums(Z, clusters, r)
  Z_stats <- cbind(Z, Z * Z)
  sums <- rowsum(Z_stats, group = clusters, reorder = FALSE)
  sums <- sums[as.character(seq_len(r)), , drop = FALSE]

  expect_equal(stats$counts, as.numeric(tabulate(clusters, nbins = r)))
  expect_equal(stats$sum_x, unname(sums[, 1L]))
  expect_equal(stats$sum_x2, unname(sums[, 2L]))
})

test_that("C++ gaussUniv cluster sums with known counts match full sums", {
  CECconfigure_runtime("fast")

  Z <- c(-2, -1, 0, 1, 3, 5)
  clusters <- c(1L, 2L, 1L, 3L, 2L, 3L)
  r <- 3L
  full <- cec_cpp_gauss_univ_cluster_sums(Z, clusters, r)
  known_counts <- cec_cpp_gauss_univ_cluster_sums_known_counts(Z, clusters, r)

  expect_equal(known_counts$sum_x, full$sum_x)
  expect_equal(known_counts$sum_x2, full$sum_x2)
})

test_that("gaussUniv cached numeric statistics match uncached optParam", {
  set.seed(104)
  Z <- c(rnorm(45, -2), rnorm(45, 2))
  phi <- phiInit(length(Z), 5)
  Z_num <- as.numeric(Z)
  Z_stats <- cbind(Z_num, Z_num * Z_num)

  uncached <- optParam_gaussUniv(Z = Z, phi = phi)
  cached <- optParam_gaussUniv(
    Z = Z,
    phi = phi,
    Z_num = Z_num,
    Z_stats = Z_stats
  )

  expect_equal(cached$nu, uncached$nu, tolerance = 1e-12)
  expect_equal(cached$m, uncached$m, tolerance = 1e-12)
  expect_equal(cached$s, uncached$s, tolerance = 1e-12)
  expect_equal(cached$varPhi, uncached$varPhi, tolerance = 1e-12)
  expect_equal(
    CECentropy_gaussUniv_from_params(cached),
    CECentropy_gaussUniv_from_params(uncached),
    tolerance = 1e-12
  )
})

test_that("gaussUniv phi matrix helper matches vectorized soft optParam", {
  set.seed(105)
  Z <- c(rnorm(35, -2), rnorm(35, 2))
  r <- 5L
  phiM <- matrix(runif(length(Z) * r), length(Z), r)
  phiM <- (1 / rowSums(phiM)) * phiM
  phi <- as.vector(phiM)
  Z_num <- as.numeric(Z)
  Z_stats <- cbind(Z_num, Z_num * Z_num)

  from_vector <- optParam_gaussUniv(
    Z = Z,
    phi = phi,
    Z_num = Z_num,
    Z_stats = Z_stats
  )
  from_matrix <- CECoptParam_gaussUniv_from_phiM(
    Z = Z,
    phiM = phiM,
    Z_stats = Z_stats
  )

  expect_equal(from_matrix$nu, from_vector$nu, tolerance = 1e-12)
  expect_equal(from_matrix$m, from_vector$m, tolerance = 1e-12)
  expect_equal(from_matrix$s, from_vector$s, tolerance = 1e-12)
  expect_equal(from_matrix$varPhi, from_vector$varPhi, tolerance = 1e-12)
  expect_equal(
    CECentropy_gaussUniv_from_params(from_matrix),
    CECentropy_gaussUniv_from_params(from_vector),
    tolerance = 1e-12
  )
})

test_that("gaussUniv optParam accepts equivalent phi matrix input", {
  set.seed(107)
  Z <- c(rnorm(40, -2), rnorm(40, 2))
  r <- 4L
  phiM <- matrix(runif(length(Z) * r), length(Z), r)
  phiM <- sweep(phiM, 1, rowSums(phiM), "/")
  phi <- as.vector(phiM)
  Z_num <- as.numeric(Z)
  Z_stats <- cbind(Z_num, Z_num * Z_num)

  from_vector <- optParam_gaussUniv(
    Z = Z,
    phi = phi,
    Z_num = Z_num,
    Z_stats = Z_stats
  )
  from_matrix <- optParam_gaussUniv(
    Z = Z,
    phi = phiM,
    Z_num = Z_num,
    Z_stats = Z_stats
  )

  expect_equal(from_matrix$nu, from_vector$nu, tolerance = 1e-12)
  expect_equal(from_matrix$m, from_vector$m, tolerance = 1e-12)
  expect_equal(from_matrix$s, from_vector$s, tolerance = 1e-12)
  expect_equal(from_matrix$varPhi, from_vector$varPhi, tolerance = 1e-12)
})

test_that("gaussUniv backend data caches invariant numeric vectors", {
  Z <- seq(-2, 2, length.out = 17)
  backend_data <- CECprepare_backend_data(Z, familyType = "gaussUniv")

  expect_equal(backend_data$Z_num, as.numeric(Z))
  expect_equal(backend_data$Z_stats[, 1], as.numeric(Z))
  expect_equal(backend_data$Z_stats[, 2], as.numeric(Z) * as.numeric(Z))
  expect_equal(length(backend_data$fingerprint_weights), length(Z))
})

test_that("integer cluster compression info matches compressed labels", {
  clusters <- c(3L, 3L, 1L, 5L, 1L)
  info <- CECcompress_integer_clusters_info(clusters, r = 5L)

  expect_equal(info$clusters, CECcompress_integer_clusters(clusters, r = 5L))
  expect_equal(info$r, max(info$clusters))
  expect_equal(info$r, 3L)
  expect_equal(info$counts, as.numeric(tabulate(info$clusters, nbins = info$r)))

  full <- c(1L, 2L, 3L, 2L)
  full_info <- CECcompress_integer_clusters_info(full, r = 3L)
  expect_equal(full_info$clusters, full)
  expect_equal(full_info$r, 3L)
  expect_equal(full_info$counts, c(1, 2, 1))
})

test_that("C++ integer cluster compression matches R reference cases", {
  clusters <- c(4L, 4L, 2L, 2L, 2L, 5L)

  CECconfigure_runtime("base")
  reference <- CECcompress_integer_clusters_info(clusters, r = 5L)
  CECconfigure_runtime("fast")
  fast <- CECcompress_integer_clusters_info(clusters, r = 5L)

  expect_equal(fast$clusters, reference$clusters)
  expect_equal(fast$r, reference$r)
  expect_equal(fast$counts, reference$counts)

  full <- c(1L, 2L, 3L, 2L, 1L)
  CECconfigure_runtime("base")
  full_reference <- CECcompress_integer_clusters_info(full, r = 3L)
  CECconfigure_runtime("fast")
  full_fast <- CECcompress_integer_clusters_info(full, r = 3L)

  expect_equal(full_fast$clusters, full_reference$clusters)
  expect_equal(full_fast$r, full_reference$r)
  expect_equal(full_fast$counts, full_reference$counts)
})

test_that("gaussUniv tagged cluster starts avoid dense phi materialization", {
  clusters <- c(2L, 1L, 2L, 3L)
  start <- CECtag_cluster_start(clusters, r = 3L)

  expect_equal(length(start), length(clusters))
  expect_true(isTRUE(attr(start, "CEC_cluster_start", exact = TRUE)))
  expect_equal(CECphi_clusters_attr(start, n = length(clusters), r = 3L), clusters)
})

test_that("gaussUniv cluster params match with provided counts", {
  set.seed(106)
  Z <- c(rnorm(35, -2), rnorm(35, 2), rnorm(20, 5))
  clusters_raw <- sample(c(1L, 3L, 5L), length(Z), replace = TRUE)
  info <- CECcompress_integer_clusters_info(clusters_raw, r = 5L)
  Z_num <- as.numeric(Z)
  Z_stats <- cbind(Z_num, Z_num * Z_num)

  without_counts <- CECoptParam_gaussUniv_from_clusters(
    Z = Z,
    clusters = info$clusters,
    compress = FALSE,
    Z_stats = Z_stats,
    r = info$r
  )
  with_counts <- CECoptParam_gaussUniv_from_clusters(
    Z = Z,
    clusters = info$clusters,
    compress = FALSE,
    Z_stats = Z_stats,
    r = info$r,
    counts = info$counts
  )

  expect_equal(with_counts$nu, without_counts$nu, tolerance = 1e-12)
  expect_equal(with_counts$m, without_counts$m, tolerance = 1e-12)
  expect_equal(with_counts$s, without_counts$s, tolerance = 1e-12)
  expect_equal(with_counts$varPhi, without_counts$varPhi, tolerance = 1e-12)
})
