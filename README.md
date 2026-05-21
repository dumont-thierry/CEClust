# CEClust

`CEClust` implements composite entropy clustering for one-dimensional,
multidimensional, categorical, and mixed data. The package is built around a
simple idea: do not trust one isolated clustering run before looking at how the
partition changes when the regularisation changes.

In practice, the two main tuning parameters are:

- `lambda`, which controls the regularisation and therefore the granularity of
  the partition;
- `C`, a density bound used by the CEC model. When `C` is poorly chosen, some
  solutions can become saturated and should not be interpreted as stable
  discoveries.

For mixed data, `C` controls only the quantitative density bound. Qualitative
variables have a separate optional floor, `Cquali`; by default `Cquali = Inf`,
so no qualitative frequency floor is applied. The reported saturation `sat`
therefore reflects quantitative density saturation, not categorical smoothing.

The recommended workflow is visual first:

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

CECplotGrid(...)         # C x lambda map
CECplotPath(...)         # REO path or 1D partition path
CECplotPartition(...)    # one selected partition
CECexplore(...)          # Shiny explorer

CECsummariseGrid(...)
CECselectStableLambdas(...)
```

By default, `CECfitBoundGrid()` keeps results in memory and does not write
checkpoint or result files. Checkpoints are an explicit option for long runs,
not part of the basic user call.

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
- `mlbench` for the breast-cancer data example;
- `mclust` for optional ICL comparisons;
- `knitr` and `rmarkdown` for vignettes.

## Principle

For fixed `C` and `lambda`, `CECclassif()` runs several random starts and keeps
the best solution according to the composite entropy criterion:

```math
H = H(\nu) + \sum_x \nu_x H(G_x \mid g_{\theta_x})
```

A single run is useful when `C` and `lambda` are already chosen. For model
selection, the package follows complete lambda trajectories and, when needed,
complete `(C, lambda)` grids. The plots are deliberately the first outputs to
look at: they show whether a candidate solution lives inside a coherent region
or is just a local accident.

## 1D Gaussian-Mixture Example

This is the main toy example. It simulates data from a Gaussian mixture, then
evaluates the same kind of `(C, lambda)` grid as the article workflow. This is
preferable to using one standard Gaussian, where artificial clusters can be
created by overfitting.

```r
library(CEClust)

set.seed(13)
n <- 100
mu <- 2 * sqrt(3) * c(-3, -1.75, -0.75, 0, 0.75, 1.75, 3)
sd_vec <- c(0.2, 1, 0.5, 1, 0.5, 1, 0.2)
pi_vec <- c(0.5, 2, 0.5, 1, 0.5, 2, 0.5)
pi_vec <- pi_vec / sum(pi_vec)

component <- sample.int(length(mu), size = n, replace = TRUE, prob = pi_vec)
Z_1d <- rnorm(n, mean = mu[component], sd = sd_vec[component])
```

Fit a full `(C, lambda)` grid on these data:

```r
partition_over_grid <- CECfitBoundGrid(
  Z = Z_1d,
  dataset_name = "gaussian_mixture_1d",
  algo_seed = 2026035,
  lambda_grid = seq(0.1, 2, by = 0.1),
  C_grid = seq(0.25, 8, by = 0.25),
  r0 = 10,
  k0 = 20,
  B = 20,
  Nshots_fresh = 50,
  Nshots_warm = 50,
  Nloop = 100,
  familyType = "gaussUniv",
  stab_algo_threshold = 0.95,
  rsi_threshold = 0.95,
  sat_threshold = 0
)
```

Start with plots. The grid shows retained and rejected cells; the path view
shows how the partition changes as `lambda` varies for a fixed `C`.

```r
CECplotGrid(partition_over_grid)

CECplotPath(
  partition_over_grid,
  C = 1,
  type = "partitions",
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.8,
  sat_threshold = 0.05,
  saturation_mark = "hatch"
)

CECplotPartition(partition_over_grid, C = 1, lambda = 1)
```

Once the data and plots are understood, fit one simple classification at
`C = 1` and `lambda = 1`:

```r
fit_1_1 <- CECclassif(
  Z = Z_1d,
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

Diagnostics should confirm what the plots suggest, not replace visual
inspection.

```r
partition_over_grid_summary <- CECsummariseGrid(partition_over_grid)
head(partition_over_grid_summary)

stable_C1 <- CECselectStableLambdas(
  partition_over_grid,
  C = 1,
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.8,
  sat_threshold = 0.05
)

stable_C1$kept_lambdas
stable_C1$summary
```

The same grid can be explored interactively in Shiny:

```r
CECexplore(partition_over_grid)
```

## Iris Without Species

Here the clustering is unsupervised: only the four quantitative iris variables
are used. The known species labels are kept only for display and interpretation.

```r
Z_iris <- iris[, 1:4]

iris_no_species <- CECfitBoundGrid(
  Z = Z_iris,
  display_data = iris,
  labels = iris$Species,
  label_name = "Species",
  dataset_name = "iris_quantitative",
  algo_seed = 2026030,
  lambda_grid = seq(0.1, 4, by = 0.5),
  C_grid = seq(1, 10, by = 2),
  r0 = 10,
  k0 = 20,
  B = 20,
  Nshots_fresh = 20,
  Nshots_warm = 20,
  Nloop = 100,
  familyType = "gaussVector"
)
```

Plot first, then inspect one selected partition or launch the Shiny explorer.

```r
CECplotGrid(iris_no_species) 
CECplotPartition(iris_no_species, shape_var = "Species",C = 1, lambda = 1)
CECexplore(iris_no_species)
```

A simple fixed-parameter classification on the same variables is:

```r
iris_fit <- CECclassif(
  Z = Z_iris,
  C = 4,
  lambda = 1,
  r0 = 10,
  Nshots = 20,
  Nloop = 100,
  familyType = "gaussVector"
)

table(iris_fit$clusters, iris$Species)
```

## Iris With Species As A Qualitative Variable

This is not the usual unsupervised iris benchmark: `Species` is now included in
the clustering variables. It is useful when you explicitly want to study a mixed
quantitative/categorical representation. With the default `Cquali = Inf`, the
qualitative variable is not treated as a saturation constraint.

```r
Z_iris_species <- data.frame(
  iris[, 1:4],
  Species = iris$Species
)

iris_with_species <- CECfitBoundGrid(
  Z = Z_iris_species,
  display_data = Z_iris_species,
  labels = iris$Species,
  label_name = "Species",
  dataset_name = "iris_with_species",
  algo_seed = 2026030,
  lambda_grid = seq(0.1, 4, by = 0.5),
  C_grid = seq(1, 10, by = 2),
  Cquali = Inf,
  r0 = 10,
  k0 = 20,
  B = 20,
  Nshots_fresh = 20,
  Nshots_warm = 20,
  Nloop = 100,
  familyType = "gaussAndDiscreteVector"
)

CECplotGrid(iris_with_species)
CECplotPartition(iris_with_species, shape_var = "Species",C = 1, lambda = 1)
CECexplore(iris_with_species)
```

The corresponding fixed-parameter classification is:

```r
iris_species_fit <- CECclassif(
  Z = Z_iris_species,
  C = 4,
  Cquali = Inf,
  lambda = 1,
  r0 = 10,
  Nshots = 20,
  Nloop = 100,
  familyType = "gaussAndDiscreteVector"
)

table(iris_species_fit$clusters, iris$Species)
```

## Uniform Example On [0, 1]

The uniform example is a stress test. There is no Gaussian-mixture structure to
recover, so the lambda path should be interpreted with caution. Here the
diagnostic thresholds are deliberately strict: keep lambda values only when
`RSI > 0.95`, `stab_algo > 0.95`, and saturation is zero.

```r
set.seed(13)
Z_uniform <- runif(1000)

uniform_grid <- CECfitBoundGrid(
  Z = Z_uniform,
  dataset_name = "uniform_0_1", 
  lambda_grid = seq(0.1, 2, by = 0.1),
  C_grid = seq(0.5, 10, by = 0.5),
  r0 = 50,
  k0 = 20,
  B = 20,
  Nshots_fresh = 20,
  Nshots_warm = 20,
  Nloop = 100,
  familyType = "gaussUniv",
  stab_algo_threshold = 0.8,
  rsi_threshold = 0.95,
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
  sat_threshold = 0,
  saturation_mark = "hatch"
)

CECplotGrid(uniform_grid)

```
Plot the phase diagram without stability or saturation criteria : 

```r
CECplotGrid(
  uniform_grid,
  stab_algo_threshold = 0 ,
  rsi_threshold = 0,
  sat_threshold = 1
)

```

The corresponding diagnostic table is:

```r
uniform_stable <- CECselectStableLambdas(
  uniform_grid,
  C = 10,
  stab_algo_threshold = 0.95,
  rsi_threshold = 0.95,
  sat_threshold = 0
)

uniform_stable$summary
```

## Further Reading

```r
vignette("minimal-api", package = "CEClust")
vignette("intro", package = "CEClust")
```
