library(CEClust)
N <- asNamespace("CEClust")
X <- matrix(c(seq(-3, -1, length.out = 12), seq(1, 3, length.out = 12)), ncol = 1)
near <- function(x, y) stopifnot(isTRUE(all.equal(x, y, tolerance = 1e-10)))
bad <- function(expr) stopifnot(inherits(try(expr, silent = TRUE), "try-error"))

# Cmin uses ML covariance and behaves correctly under affine rescaling.
g <- CECinitialGrid(X)
near(g$Cmin, 1 / sqrt(2 * pi * mean((X - mean(X))^2)))
near(g$C_grid, g$Cmin * 2^(1:10))
near(CECinitialGrid(3 * X + 11)$Cmin, g$Cmin / 3)
Y <- cbind(X, sin(X))
near(CECinitialGrid(Y)$Cmin, 1 / (2 * pi * sqrt(det(cov(Y) * 23 / 24))))
bad(CECinitialGrid(cbind(X, X)))
bad(CECinitialGrid(X, c(1, 3, 4)))
bad(CECdiagramControl("adaptive", typo = 1))
bad(CECdiagramControl("bands", typo = 1))
bad(CECdiagram(X, "bands", C_grid = c(1, 2)))
stopifnot(CECdiagramControl("bands")$max_anchors == 10,
  CECdiagramControl("bands")$max_projection_queries == 200,
  !CECdiagramControl("adaptive")$adaptive_control$extend_left,
  N$CECadaptive_control(list(), 0)$extend_left)

boots <- list(rep(seq_len(12), each = 2), rep(13:24, each = 2))
ac <- CECdiagramControl("adaptive", optimizer_control = list(Rmax = 2, Ninit = 1,
  max_iter = 8, max_batches = 1), adaptive_control = list(max_cells = 60,
  max_rounds = 2, max_extensions = 0, diagnostic_strategy = "cells"))
for (penalty in c("none", "bic")) {
  path <- tempfile("unified-grid-")
  fit <- CECdiagram(X, "adaptive", L = 2, B = 2, seed = 19, penalty = penalty,
    C_grid = c(1, 5), lambda_grid = c(.6, 1.4), bootstrap_indices = boots,
    control = ac, directory = path)
  stopifnot(inherits(fit, "CEC_adaptive_grid"), fit$meta$penalty == penalty,
    fit$adaptive$extensions['left'] == 0, file.exists(file.path(path, "diagram.rds")))
  direct <- do.call(N$CECfit_adaptive_grid, c(list(Z = X, k0 = 2, B = 2,
    algo_seed = 19, penalty = penalty, C_grid = c(1, 5), lambda_grid = c(.6, 1.4),
    bootstrap_indices = boots), ac))
  near(fit$summary, direct$summary)
  again <- CECdiagram(X, "adaptive", L = 2, B = 2, seed = 19, penalty = penalty,
    C_grid = c(1, 5), lambda_grid = c(.6, 1.4), bootstrap_indices = boots,
    control = ac, directory = path, resume = TRUE)
  near(fit$summary, again$summary)
  bad(CECdiagram(X, "adaptive", control = ac, directory = path))
}

# Automatic grid and fixed-left boundary with extension budget available.
auto <- CECdiagram(X, L = 2, B = 1, lambda_grid = c(10, 20), control = list(
  optimizer_control = list(Rmax = 1, rinit_max = 1, Ninit = 1, max_iter = 2, max_batches = 1),
  adaptive_control = list(max_cells = 80, max_extensions = 1, max_rounds = 2)))
near(auto$adaptive$initial_C, g$C_grid)
near(attr(auto, "CECdiagram")$initialization$Cmin, g$Cmin)
stopifnot(auto$adaptive$extensions['left'] == 0)

# Complete bands, BIC, checkpoint and worker routing through the unified namespace.
bc <- CECdiagramControl("bands", orders = 1:2, starts = 1, extra_starts = 0,
  max_order = 4, max_iter = 12, max_anchors = 2, C_points = 3, max_C_points = 7,
  refine_rounds = 0, projection_rounds = 0, max_projection_queries = 4)
for (penalty in c("none", "bic")) {
  path <- tempfile("unified-bands-")
  a <- CECdiagram(X, "bands", L = 2, B = 2, workers = 1, seed = 19,
    bootstrap_indices = boots, penalty = penalty, control = bc, directory = path)
  b <- CECdiagram(X, "bands", L = 2, B = 2, workers = 2, seed = 19,
    bootstrap_indices = boots, penalty = penalty, control = bc,
    directory = tempfile("unified-parallel-"))
  stopifnot(bands_validate(a)$passed, bands_validate(b)$passed)
  near(a$tiles, b$tiles)
  for (name in names(a$runs)) near(names(a$runs[[name]]$candidates), names(b$runs[[name]]$candidates))
  r <- CECdiagram(X, "bands", L = 2, B = 2, workers = 1, seed = 19,
    bootstrap_indices = boots, penalty = penalty, control = bc, directory = path, resume = TRUE)
  near(a$tiles, r$tiles)
  bad(CECdiagram(X + 1, "bands", L = 2, B = 2, seed = 19,
    bootstrap_indices = boots, penalty = penalty, control = bc, directory = path, resume = TRUE))
}
cat("Unified interface, Cmin, both penalties, resumption and parallel bands: passed\n")
