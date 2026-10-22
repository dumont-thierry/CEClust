% Maintained by hand; see CONTRIBUTING.md.
\name{CECextractPartition}
\alias{CECextractPartition}
\title{Extract the labels at a grid cell}
\description{
Obtain the representative partition stored at one cell. Cells need not pass the diagnostic filter to be inspected.
}
\usage{
CECextractPartition(x, C, lambda)
}
\arguments{
\item{x}{A grid returned by \code{CECfitBoundGrid()} or \code{CECfitPreset()}. For a consolidated one-C diagnostic, use its \code{grid} component.}
\item{C, lambda}{Grid coordinates. The nearest stored C and then the nearest stored lambda are used; no new optimization or interpolation is performed. Use coordinates from \code{CECsummariseGrid(x)} to select an exact cell.}
}
\value{
An integer vector with one cluster label per original observation, in the original row order.
}
\details{
Labels have no intrinsic meaning and may be renumbered in a plot.
The partition is the representative selected by the fitting engine, not an
individual bootstrap output. Inspect \code{keep}, \code{measured} and
\code{sat} in the summary before interpreting it as retained. Use
\code{CECstability()} to compare partitions up to relabeling.
}
\seealso{\link{CECsummariseGrid}, \link{CECextractBestPartitions}, \link{CECplotPartition}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
g <- CECfitBoundGrid(Z = X, C_grid = c(0.2, 1),
  lambda_grid = c(0.6, 1, 1.5), optimization = "em_bank",
  penalty = "none", k0 = 3, B = 3, algo_seed = 14,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  n_cores = 1)
z <- CECextractPartition(g, C = 1, lambda = 1)
table(z)
}
