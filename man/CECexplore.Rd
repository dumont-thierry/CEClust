% Maintained by hand; see CONTRIBUTING.md.
\name{CECexplore}
\alias{CECexplore}
\title{Open the interactive grid explorer}
\description{
Inspect an already fitted grid in a Shiny application. Select cells, view their partitions and vary display filters without refitting the clustering.
}
\usage{
CECexplore(x, ...)
}
\arguments{
\item{x}{A grid returned by \code{CECfitBoundGrid()} or \code{CECfitPreset()}. For a consolidated one-C diagnostic, use its \code{grid} component.}
\item{...}{\code{launch.browser=TRUE} opens a browser; set FALSE to launch without opening one.}
}
\value{
The return value of \code{shiny::runApp()}; the call normally blocks until the application is stopped.
}
\details{
Requires the optional package \pkg{shiny}. The application operates on
the stored grid, partitions and diagnostic values. Changing thresholds cannot
compute an RSI that was skipped during fitting. To enlarge the measurement
policy for an EM grid, use \code{CECrescoreGrid()} before reopening it.
Stop the application in the R session to return to the console.
}
\seealso{\link{CECplotGrid}, \link{CECrescoreGrid}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
g <- CECfitBoundGrid(Z = X, C_grid = c(0.2, 1),
  lambda_grid = c(0.6, 1, 1.5), optimization = "em_bank",
  penalty = "none", k0 = 3, B = 3, algo_seed = 14,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  n_cores = 1)
if (interactive() && requireNamespace("shiny", quietly = TRUE)) {
  CECexplore(g, launch.browser = TRUE)
}
}
