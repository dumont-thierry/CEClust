library(CEClust)

make_sigma_2d <- function(sd_major, sd_minor = sd_major, angle = 0) {
  rot <- matrix(
    c(cos(angle), -sin(angle), sin(angle), cos(angle)),
    nrow = 2,
    byrow = TRUE
  )
  rot %*% diag(c(sd_major^2, sd_minor^2), nrow = 2) %*% t(rot)
}

simulate_mixture_2d <- function(n, seed = 13) {
  set.seed(seed)

  weights <- c(0.45, 0.25, 0.30)
  means <- list(c(-2, 2), c(0, 0), c(5, -6))
  covariances <- list(
    0.16 * diag(2),
    5.76 * diag(2),
    matrix(c(8.5, 7.5, 7.5, 8.5), nrow = 2)
  )

  z_true <- sample.int(3, size = n, replace = TRUE, prob = weights)
  Z <- matrix(NA_real_, nrow = n, ncol = 2)

  for (j in seq_along(weights)) {
    idx <- which(z_true == j)
    if (!length(idx)) {
      next
    }

    eig <- eigen(covariances[[j]], symmetric = TRUE)
    transform <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)), nrow = 2)
    noise <- matrix(rnorm(2 * length(idx)), ncol = 2)
    Z[idx, ] <- matrix(means[[j]], nrow = length(idx), ncol = 2, byrow = TRUE) +
      noise %*% t(transform)
  }

  colnames(Z) <- c("z1", "z2")
  list(
    Z = as.data.frame(Z),
    z_true = factor(z_true, labels = c("compact", "central", "elongated"))
  )
}

mixture_2d <- simulate_mixture_2d(n = 300, seed = 13)

plot(
  mixture_2d$Z$z1,
  mixture_2d$Z$z2,
  col = as.integer(mixture_2d$z_true),
  pch = 19,
  xlab = expression(z[1]),
  ylab = expression(z[2])
)

grid_2d <- CECfitBoundGrid(
  Z = mixture_2d$Z,
  display_data = data.frame(mixture_2d$Z, true_component = mixture_2d$z_true),
  labels = mixture_2d$z_true,
  label_name = "Generating component",
  dataset_name = "synthetic_mixture_2d",
  algo_seed = 202606,
  lambda_grid = seq(3 / 20, 3, length.out = 12),
  C_grid = seq(1 / 10, 2, length.out = 12),
  r0 = 10,
  k0 = 20,
  B = 10,
  Nshots_fresh = 10,
  Nshots_warm = 10,
  Nloop = 100,
  familyType = "gaussVector",
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.95,
  sat_threshold = 0,
  n_cores = 2,
  compact_results = TRUE
)

CECplotGrid(
  grid_2d,
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.95,
  sat_threshold = 0
)

CECplotPartition(grid_2d, C = grid_2d$C_grid[8], lambda = grid_2d$lambda_grid[6])

summary_2d <- CECsummariseGrid(grid_2d)
head(summary_2d)

if (interactive()) {
  CECexplore(grid_2d)
}
