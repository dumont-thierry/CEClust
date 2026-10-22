% Maintained by hand; see CONTRIBUTING.md.
\name{CECplotPartition}
\alias{CECplotPartition}
\title{Inspect the partition at one cell}
\description{
Plot a selected partition against one-dimensional data or two displayed coordinates. The fitting coordinates and the display coordinates can differ.
}
\usage{
CECplotPartition(x, C = NULL, lambda = NULL, x_var = NULL, y_var = NULL, shape_var = "None", 
    ...)
}
\arguments{
\item{x}{A grid returned by \code{CECfitBoundGrid()} or \code{CECfitPreset()}. For a consolidated one-C diagnostic, use its \code{grid} component.}
\item{C, lambda}{Grid coordinates. The nearest stored C and then the nearest stored lambda are used; no new optimization or interpolation is performed. Use coordinates from \code{CECsummariseGrid(x)} to select an exact cell.}
\item{x_var, y_var}{Names of displayed quantitative variables, or \code{"PCA1"}, \code{"PCA2"}, etc. Defaults use the first available display choices.}
\item{shape_var}{Name of a qualitative display variable used for point shapes, or \code{"None"}. It need not be a fitting feature.}
\item{...}{For multivariate displays: \code{cluster_selection="All"} or a group selection, and \code{point_cex=1.05}. For one-dimensional displays: \code{hist_breaks}, \code{histogram_col}, \code{density_col}, \code{main}.}
}
\value{
Plotting information returned invisibly; NULL if the selected data cannot be displayed.
}
\details{
If coordinates are omitted, the display selects a retained cell with smallest REO, taking the first tie. If no cell is retained, it uses a finite-REO cell instead. This is a display default, not a rule for choosing a final clustering.

Coordinates select the nearest stored cell. PCA here is a display transformation,
not an instruction to refit on two axes. Display group labels can be reordered
by weight; use \code{CECextractPartition()} for the stored labels. A plot can
show a rejected cell for comparison, so always check its diagnostics.
}
\seealso{\link{CECextractPartition}, \link{CECplotGrid}, \link{CECexplore}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
g <- CECfitBoundGrid(Z = X, C_grid = c(0.2, 1),
  lambda_grid = c(0.6, 1, 1.5), optimization = "em_bank",
  penalty = "none", k0 = 3, B = 3, algo_seed = 14,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  n_cores = 1)
CECplotPartition(g, C = 1, lambda = 1, x_var = "x", y_var = "y")
}
