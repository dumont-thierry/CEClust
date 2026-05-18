# CEClust

`CEClust` implements composite entropy clustering for Gaussian, categorical,
and mixed data. The package is designed for users who want both:

- a fixed-`lambda` clustering routine with repeated random starts;
- a diagnostic workflow that studies how the selected partition evolves over a
  grid of regularisation values.

The package includes helpers for:

- selecting a conservative runtime configuration;
- fitting linked forward/backward lambda paths;
- fitting complete `C x lambda` grids;
- summarising bootstrap stability;
- extracting one representative partition per lambda;
- detecting or repairing changes along a lambda trajectory;
- visualising static outputs and launching the Shiny explorer.

The minimal public workflow is:

```r
CECclassif(...)          # one lambda
CECfitLambdaGrid(...)    # one C, several lambdas
CECfitBoundGrid(...)     # several C values and several lambdas
CECfitPreset("iris")     # built-in demonstration workflows

CECselectStableLambdas(...)
CECsummariseGrid(...)

CECplotGrid(...)
CECplotPartition(...)
CECplotPath(...)
CECexplore(...)          # Shiny explorer
```

## Package goal

The main scientific goal of `CEClust` is not only to produce one clustering
solution, but also to help the user understand how stable that solution is when
the regularisation parameter changes. In practice, a typical analysis proceeds
in two stages:

1. fit repeated-shot models for a sequence of `lambda` values;
2. identify stable regions where the recovered partition is robust.

This makes the package useful both for exploratory analyses and for more
structured methodological studies.

## Installation

From GitHub:

```r
remotes::install_github("dumont-thierry/CEClust")
```

From a local source folder:

```r
install.packages(".", repos = NULL, type = "source")
```

## Dependencies

Core dependencies are installed automatically with the package, including
`mvtnorm`, `clue`, `Rcpp`, and base R recommended packages used for runtime and
graphics support.

Optional packages extend the workflow:

- `ggplot2` for PCA-based graphics and example data;
- `shiny` for the interactive explorer;
- `mlbench` for the breast-cancer preset;
- `mclust` for optional ICL comparisons;
- `knitr` and `rmarkdown` for the vignette and README rendering;
- `pkgload` for development-mode parallel workers when the package has not yet
  been installed.

## Method overview

For a fixed `lambda`, `CECclassif()` runs several initialisations and keeps the
best solution according to the composite entropy criterion:

```text
H = H(nu) + sum_x nu_x H(G_x | g_theta_x)
```

where:

- `H(nu)` measures the entropy of cluster weights;
- the conditional term measures within-cluster model fit;
- `lambda` controls the balance between these components.

The linked-lambda workflow extends this by following the fitted partition
forward and backward along a grid of `lambda` values. Bootstrap paths are then
projected back to the original data to quantify stability.

## Minimal data example

The package provides a small simulator for mixed data:

```r
library(CEClust)

Z <- simulate_multidim_benchmark_data(
  n = 120,
  p_num = 3,
  p_fac = 2,
  seed = 1
)

str(Z)
```

For a ready-made mixed example based on a real data set:

```r
Z_diamonds <- CECsample_diamonds(n_max = 100, seed = 1)
str(Z_diamonds)
```

## Basic workflow

### 1. Configure the runtime

```r
runtime <- CECconfigure_runtime("auto", max_cores = 4)
runtime
```

Use `"base"` if you want a conservative single-core pure-R setup. The returned
`n_cores` value can be passed directly to the lambda-grid diagnostic workflow.

### 2. Fit a model for one lambda

```r
fit <- CECclassif(
  Z = Z,
  lambda = 0.8,
  C = 10,
  r0 = 6,
  Nshots = 5,
  Nloop = 30,
  familyType = "gaussAndDiscreteVector"
)

fit$REO
table(fit$clusters)
```

This is the simplest entry point when you already have a candidate `lambda`
value in mind.

### 3. Diagnose a lambda grid

```r
lambda_diag <- CECfitLambdaGrid(
  Z = Z,
  lambda_grid = seq(0.2, 1.4, by = 0.2),
  k0 = 4,
  B = 4,
  C = 10,
  r0 = 6,
  Nshots_fresh = 2,
  Nshots_warm = 1,
  Nloop = 20,
  familyType = "gaussAndDiscreteVector",
  seed = 1,
  silent = TRUE,
  verbose = FALSE,
  n_cores = runtime$n_cores,
  checkpoint_dir = FALSE,
  auto_checkpoint = FALSE,
  show_progress = FALSE
)

head(lambda_diag$summary)
```

The `summary` component contains, for each `lambda`, the fixed-sample stability,
bootstrap stability, criterion summaries, and information about the best
initialisation origin.

### 4. Smooth and identify stable regions

```r
lambda_diag <- CECaddSmoothedDiagnostics(lambda_diag, k = 3)

stable <- CECidentifyStableLambdas(
  lambda_diag,
  rule = "ratio",
  ratio_threshold = 0.8,
  min_consecutive = 2
)

stable$lambda_min
stable$summary
```

`stable$lambda_min` gives the first lambda entering a stable region under the
chosen rule. The full summary is often more informative than a single selected
value because it shows whether stability persists across an interval.

If you also want to exclude lambdas for which too much of the sample lives in
clusters saturating the density bound `C`, extract one representative partition
per lambda and feed it back into the stability filter:

```r
best_parts <- CECextractBestPartitions(
  lambda_diag,
  source = "fixed",
  criterion = "projected_H",
  Z = Z
)

stable_sat <- CECidentifyStableLambdas(
  lambda_diag,
  rule = "ratio",
  ratio_threshold = 0.8,
  stabB_threshold = 0.8,
  sat_threshold = 0.8,
  best_partitions_obj = best_parts,
  min_consecutive = 1
)

stable_sat$summary
```

### 5. Extract one partition per lambda

```r
best_parts <- CECextractBestPartitions(
  lambda_diag,
  source = "fixed",
  criterion = "projected_H",
  Z = Z
)

head(best_parts$summary)
```

This step converts the detailed diagnostic object into a compact trajectory of
representative partitions.

### 6. Check coherence or detect changes

```r
coherence <- CECcheckBestPartitionCoherence(best_parts, Z)
changes <- CECdetectPartitionChanges(best_parts)

coherence$all_coherent
changes$change_lambdas
```

If needed, you can try a repair pass:

```r
best_parts_repaired <- CECrepairBestPartitionTrajectory(best_parts, Z)
best_parts_repaired$n_repairs
```

## Interpreting the main outputs

Some components are especially useful in practice:

- `lambda_diag$summary`: per-lambda stability and criterion diagnostics;
- `stable$summary`: stability rule applied to each lambda;
- `best_parts$summary`: representative partition chosen at each lambda;
- `coherence$summary`: local coherence checks between neighbouring lambdas;
- `changes$summary`: where the partition effectively changes along the path.

In many analyses, the most informative workflow is:

1. inspect `lambda_diag$summary`;
2. derive stable regions with `CECidentifyStableLambdas()`;
3. extract partitions with `CECextractBestPartitions()`;
4. inspect coherence and change points.

## Visualisation tools

Static summaries:

```r
plotCECdiagnoseLambdaGrid(lambda_diag)
```

One-dimensional inspection for scalar numeric data:

```r
Z_1d <- rnorm(80)

diag_1d <- CECfitLambdaGrid(
  Z = Z_1d,
  lambda_grid = seq(0.2, 1.0, by = 0.2),
  k0 = 3,
  B = 2,
  C = 10,
  r0 = 5,
  Nshots_fresh = 1,
  Nshots_warm = 1,
  Nloop = 20,
  familyType = "gaussUniv",
  seed = 1,
  silent = TRUE,
  verbose = FALSE,
  n_cores = 1L,
  checkpoint_dir = FALSE,
  auto_checkpoint = FALSE,
  show_progress = FALSE
)

best_1d <- CECextractBestPartitions(diag_1d, source = "fixed", Z = Z_1d)
plotCECBestPartitions1D(Z_1d, best_1d)
```

PCA-based trajectory display for multivariate numeric data:

```r
diag_iris <- CECfitLambdaGrid(
  iris[, -5],
  lambda_grid = seq(0.01, 0.8, by = 0.1),
  k0 = 3,
  B = 2,
  C = 10,
  r0 = 4,
  Nshots_fresh = 1,
  Nshots_warm = 1,
  Nloop = 20,
  familyType = "gaussVector",
  seed = 1,
  silent = TRUE,
  verbose = FALSE,
  n_cores = 1L,
  checkpoint_dir = FALSE,
  auto_checkpoint = FALSE,
  show_progress = FALSE
)

best_iris <- CECextractBestPartitions(diag_iris, source = "fixed", Z = iris[, -5])
path_df <- build_partition_path_df(iris[, -5], best_iris)
plot_partition_path_pca(path_df)
```

Interactive explorer:

```r
plotCECdiagnoseInteractive(best_parts, Z)
```

## Remarks and limitations

- `build_partition_path_df()` requires at least two numeric variables because it
  relies on a two-dimensional PCA representation.
- The interactive explorer requires `shiny` and `ggplot2`.
- Multi-core runs are most useful for large lambda grids or repeated bootstrap
  analyses; for quick checks, the pure-R single-core mode is easier to debug.
- The choice of `lambda` remains a modelling decision. The package helps assess
  stability, but it does not replace methodological judgement.

## Further reading

For a slightly longer walkthrough with the same core workflow, see:

```r
vignette("intro", package = "CEClust")
```
