% Maintained by hand; see CONTRIBUTING.md.
\name{CECcompactGrid}
\alias{CECcompactGrid}
\title{Remove redundant intermediate grid data}
\description{
Compact a grid while preserving its representative partitions, summary and visual inspection workflow.
}
\usage{
CECcompactGrid(x, keep_lambda_details = FALSE)
}
\arguments{
\item{x}{A grid returned by \code{CECfitBoundGrid()} or \code{CECfitPreset()}. For a consolidated one-C diagnostic, use its \code{grid} component.}
\item{keep_lambda_details}{Whether to preserve full historical per-lambda diagnostic details. Default FALSE.}
}
\value{
A compact grid with its classes retained.
}
\details{
Dense membership matrices and historical per-task details are removed.
Consolidated grids are already compact in their display components; their
\code{banks}, \code{parts}, \code{winners} and bootstrap indices remain, so
\code{CECrescoreGrid()} remains available for EM banks. This function is not a
request to discard banks or to anonymize the stored observations. Save the
returned object rather than expecting the input to change.
}
\seealso{\link{CECresults}, \link{CECrescoreGrid}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
g <- CECfitBoundGrid(Z = X, C_grid = c(0.2, 1),
  lambda_grid = c(0.6, 1, 1.5), optimization = "em_bank",
  penalty = "none", k0 = 3, B = 3, algo_seed = 14,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  n_cores = 1)
small <- CECcompactGrid(g)
stopifnot(identical(CECextractPartition(g, 1, 1),
                    CECextractPartition(small, 1, 1)))
}
