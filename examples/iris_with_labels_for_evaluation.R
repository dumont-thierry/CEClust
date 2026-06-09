library(CEClust)

Z_iris <- iris[, 1:4]

iris_grid_labeled <- CECfitBoundGrid(
  Z = Z_iris,
  display_data = iris,
  labels = iris$Species,
  label_name = "Species",
  dataset_name = "iris_labels_for_evaluation",
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

CECplotGrid(iris_grid_labeled)
CECplotPartition(iris_grid_labeled, shape_var = "Species")

iris_labeled_summary <- CECsummariseGrid(iris_grid_labeled)
head(iris_labeled_summary)

if (interactive()) {
  CECexplore(iris_grid_labeled)
}
