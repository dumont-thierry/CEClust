library(CEClust)

Z_iris <- iris[, 1:4]

iris_grid <- CECfitBoundGrid(
  Z = Z_iris,
  dataset_name = "iris_without_labels",
  algo_seed = 202606,
  lambda_grid = seq(0.4, 2.4, by = 0.2),
  C_grid = seq(1, 5, by = 0.5),
  r0 = 10,
  k0 = 20,
  B = 10,
  Nshots_fresh = 10,
  Nshots_warm = 10,
  Nloop = 100,
  familyType = "gaussVector",
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.9,
  sat_threshold = 0.05,
  n_cores = 2,
  compact_results = TRUE
)

CECplotGrid(iris_grid)
CECplotPartition(iris_grid)

iris_summary <- CECsummariseGrid(iris_grid)
head(iris_summary)

if (interactive()) {
  CECexplore(iris_grid)
}
