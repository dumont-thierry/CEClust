% Maintained by hand; see CONTRIBUTING.md.
\name{CECplotPath}
\alias{CECplotPath}
\title{Plot a diagnostic slice or one-dimensional partition path}
\description{
Inspect the group-count path at one C, or the evolution of one-dimensional allocations with rejection bands.
}
\usage{
CECplotPath(x, C = NULL, lambda = NULL, type = c("summary", "partitions"), lambda_stable_obj = NULL, 
    stab_algo_threshold = 0.8, rsi_threshold = 0.8, sat_threshold = 0.1, use_smoothed = FALSE, 
    min_consecutive = 1, ...)
}
\arguments{
\item{x}{A grid returned by \code{CECfitBoundGrid()} or \code{CECfitPreset()}. For a consolidated one-C diagnostic, use its \code{grid} component.}
\item{C}{C slice, resolved to the nearest stored value. If omitted, the display-default cell supplies it.}
\item{lambda}{Highlighted lambda for type="summary". Ignored for type="partitions".}
\item{type}{\code{"summary"} for a group-count slice; \code{"partitions"} for an allocation map, supported only for one-dimensional numeric fitting data.}
\item{lambda_stable_obj}{Optional output of \code{CECselectStableLambdas()} for the same C. Supplying it makes the path use exactly that filter.}
\item{stab_algo_threshold, rsi_threshold, sat_threshold}{Historical defaults 0.8, 0.8, 0.1 used if type="partitions" computes its own selector. These defaults do not inherit the grid thresholds; prefer passing lambda_stable_obj explicitly.}
\item{use_smoothed, min_consecutive}{Options for the computed selector. New solver grids require use_smoothed=FALSE.}
\item{...}{For a summary path, \code{main}. For the one-dimensional allocation plot, arguments of the historical partition-path renderer, such as its axis/title controls; use defaults unless reproducing an established display.}
}
\value{
The plotted summary or partition-path data, invisibly.
}
\details{
The summary path can be drawn for multivariate data, since it only shows
the group count. An allocation path requires one-dimensional fitting data.
Grey regions and saturation hatching reflect the supplied selector; they are
not additional statistical tests. If coordinates are omitted, the display selects a retained cell with smallest REO, taking the first tie. If no cell is retained, it uses a finite-REO cell instead. This is a display default, not a rule for choosing a final clustering.
}
\seealso{\link{CECselectStableLambdas}, \link{CECplotGrid}, \link{CECplotPartition}.}
\examples{
set.seed(4)
z <- c(rnorm(25, -2, 0.5), rnorm(25, 2, 0.5))
g <- CECfitBoundGrid(Z = z, C_grid = 1, lambda_grid = c(0.6, 1, 1.5),
  optimization = "em_bank", k0 = 3, B = 2, algo_seed = 9,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 20))
CECplotPath(g, C = 1, type = "summary")
keep <- CECselectStableLambdas(g, C = 1)
CECplotPath(g, C = 1, type = "partitions", lambda_stable_obj = keep)
}
