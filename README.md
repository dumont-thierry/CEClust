# CEClust

`CEClust` implements composite entropy clustering unidimensional, multidimentional, categorical,
and mixed data. The package is built around a simple idea: do not trust one
isolated clustering run before looking at how the partition changes when the
regularisation changes.

In practice, the two tuning parameters are:

- `lambda`, which controls the regularisation and therefore the granularity of
  the partition;
- `C`, a density bound used by the CEC model. When `C` is poorly chosen, some
  solutions can become saturated and should not be interpreted as stable
  discoveries.

The recommended workflow is therefore visual first:

1. fit a grid of candidate values `(C, lambda)`;
2. plot the grid and the lambda path before reading detailed diagnostics;
3. inspect one selected partition;
4. open the Shiny app when the grid is large or multivariate;
5. only then summarise stability, RSI, saturation, and selected lambdas.

The minimal public API is:

```r
CECclassif(...)          # one C and one lambda
CECfitLambdaGrid(...)    # one C and several lambdas
CECfitBoundGrid(...)     # several C values and several lambdas
CECfitPreset("iris")     # built-in workflows

CECplotGrid(...)         # C x lambda map
CECplotPath(...)         # REO path or 1D partition path
CECplotPartition(...)    # one selected partition
CECexplore(...)          # Shiny explorer

CECsummariseGrid(...)
CECselectStableLambdas(...)
```

By default, `CECfitBoundGrid()` and `CECfitPreset()` keep results in memory and
do not write checkpoint or result files. Checkpoints are an explicit option for
long runs, not part of the basic user call.

## Installation

Current development branch:

```r
remotes::install_github(
  "dumont-thierry/CEClust",
  ref = "minimal-api",
  build_vignettes = TRUE
)
```
 
From a local source folder:

```r
install.packages(".", repos = NULL, type = "source")
```

Optional packages extend the workflow:

- `shiny` for `CECexplore()`;
- `ggplot2` for some multivariate displays;
- `mlbench` for the breast-cancer preset;
- `mclust` for optional ICL comparisons;
- `knitr` and `rmarkdown` for vignettes.

## Principle

For fixed `C` and `lambda`, `CECclassif()` runs several random starts and keeps
the best solution according to the composite entropy criterion:

\[
H = H(\nu) + \sum_x \nu_x H(G_x \mid g_{\theta_x})
\]

A single run is useful when `C` and `lambda` are already chosen. For model
selection, the package follows complete lambda trajectories and, when needed,
complete `(C, lambda)` grids. The plots are deliberately the first outputs to
look at: they show whether a candidate solution lives inside a coherent region
or is just a local accident.

## 1D Gaussian-mixture example

This is the main toy example. It simulates data from the package's Gaussian
mixture model, then evaluates the same kind of grid as the article workflow.
This replaces the misleading example where a single standard Gaussian could be
split into several artificial clusters.

```r
library(CEClust)

article_grid <- CECfitPreset(
  "gaussian_1d",
  n = 100,
  data_seed = 13,
  algo_seed = 2026035,
  lambda_grid = seq(0.1, 1.8, by = 0.1),
  C_grid = c(0.5, 1, 2, 3, 3.5, 4, 4.5, 5, 6, 8),
  r0 = 10,
  k0 = 20,
  B = 20,
  Nshots_fresh = 10,
  Nshots_warm = 10,
  Nloop = 100,
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.8,
  sat_threshold = 0.05
)
```

### Start with plots

The first plot is the `(C, lambda)` grid. Kept cells are coloured by REO;
rejected cells are hatched according to instability, RSI, or saturation.

```r
CECplotGrid(article_grid)
```

For one-dimensional data, the lambda path can be read as a sequence of
partitions. Horizontal bands mark values of `lambda`; grey or hatched regions
are values rejected by the diagnostic rules.

```r
CECplotPath(
  article_grid,
  C = 1,
  type = "partitions",
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.8,
  sat_threshold = 0.05,
  saturation_mark = "hatch"
)
```

Then inspect the selected partition itself.

```r
CECplotPartition(article_grid, C = 1, lambda = 1)
```

### Same data, one simple classification

Once the data and plots are understood, a single classification at `C = 1` and
`lambda = 1` is straightforward:

```r
fit_1_1 <- CECclassif(
  Z = article_grid$Z,
  C = 1,
  lambda = 1,
  r0 = 10,
  Nshots = 50,
  Nloop = 100,
  familyType = "gaussUniv"
)

table(fit_1_1$clusters)
fit_1_1$REO
```

### Then summarise diagnostics

Diagnostics should confirm what the plots suggest, not replace visual
inspection.

```r
article_summary <- CECsummariseGrid(article_grid)
head(article_summary)

stable_C1 <- CECselectStableLambdas(
  article_grid,
  C = 1,
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.8,
  sat_threshold = 0.05
)

stable_C1$kept_lambdas
stable_C1$summary
```

### Open the Shiny grid explorer

The same grid can be explored interactively. This is the easiest way to move
between the grid, selected partitions, variables, clusters, and lambda paths.

```r
CECexplore(article_grid)
```

## Iris without and with Species

The iris preset keeps the same API but changes the data type. First run the
classification on the four quantitative variables only.

```r
iris_no_species <- CECfitPreset(
  "iris",
  include_species_as_feature = FALSE,
  quantitative_representation = "raw",
  lambda_grid = seq(0.1, 4, by = 0.5),
  C_grid = seq(1, 10, by = 2),
  r0 = 10,
  k0 = 20,
  B = 20,
  Nshots_fresh = 20,
  Nshots_warm = 20,
  Nloop = 100
)

CECplotGrid(iris_no_species)
CECplotPartition(iris_no_species)
CECexplore(iris_no_species)
```

Then include `Species` as a qualitative variable in the clustering data. This
is not the usual unsupervised iris benchmark; it is useful when you explicitly
want to study a mixed quantitative/categorical representation.

```r
iris_with_species <- CECfitPreset(
  "iris",
  include_species_as_feature = TRUE,
  quantitative_representation = "raw",
  lambda_grid = seq(0.1, 4, by = 0.5),
  C_grid = seq(1, 10, by = 2),
  r0 = 10,
  k0 = 20,
  B = 20,
  Nshots_fresh = 20,
  Nshots_warm = 20,
  Nloop = 100
)

CECplotGrid(iris_with_species)
CECplotPartition(iris_with_species, shape_var = "Species")
CECexplore(iris_with_species)
```

## Uniform example on [0, 1]

The uniform example is a stress test. There is no Gaussian-mixture structure to
recover, so the lambda path should be interpreted with caution. Here the
diagnostic thresholds are deliberately strict: keep lambda values only when
`RSI > 0.9`, `stab_algo > 0.8`, and saturation is zero.

```r
set.seed(13)
Z_uniform <- runif(1000)

uniform_grid <- CECfitBoundGrid(
  Z = Z_uniform,
  dataset_name = "uniform_0_1",
  lambda_grid = seq(0.05, 2, by = 0.05),
  C_grid = 10,
  r0 = 50,
  k0 = 20,
  B = 20,
  Nshots_fresh = 20,
  Nshots_warm = 20,
  Nloop = 100,
  familyType = "gaussUniv",
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.9,
  sat_threshold = 0
)
```

Plot the partition for each lambda value at `C = 10`. Hatching marks lambdas
that fail the requested rules, including the `sat = 0` rule.

```r
CECplotPath(
  uniform_grid,
  C = 10,
  type = "partitions",
  stab_algo_threshold = 0.95,
  rsi_threshold = 0.95,
  sat_threshold = 1.1,
  saturation_mark = "hatch"
)

CECplotGrid(uniform_grid)
```

The corresponding diagnostic table is:

```r
uniform_stable <- CECselectStableLambdas(
  uniform_grid,
  C = 10,
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.9,
  sat_threshold = 0
)

uniform_stable$summary
```

## Presets

Presets are named entry points for examples that should remain easy to launch:

```r
CECfitPreset("gaussian_1d")
CECfitPreset("iris")
CECfitPreset("breast_cancer")
```

They are intentionally thin wrappers around `CECfitBoundGrid()`. Arguments such
as `lambda_grid`, `C_grid`, `k0`, `B`, `Nshots_fresh`, `Nloop`, `n_cores`, and
checkpoint options can still be passed explicitly.

## Further reading

```r
vignette("minimal-api", package = "CEClust")
vignette("intro", package = "CEClust")
```
