% Maintained by hand; see CONTRIBUTING.md.
\name{CECaddSmoothedDiagnostics}
\alias{CECaddSmoothedDiagnostics}
\title{Add moving averages to historical lambda diagnostics}
\description{
Compatibility helper for descriptive smoothing of an existing lambda-summary table. It does not refit the algorithm or create new independent comparisons.
}
\usage{
CECaddSmoothedDiagnostics(lambda_diag, k = 3, ratio_col = "stab_ratio_mean", stab0_col = "stab0_mean", 
    stabB_col = "stabB_mean", H_col = "meanH_B")
}
\arguments{
\item{lambda_diag}{Historical lambda diagnostics containing the requested summary columns.}
\item{k}{Moving-average window in grid entries, default 3.}
\item{ratio_col, stab0_col, stabB_col, H_col}{Names of the ratio, algorithmic stability, bootstrap stability and score columns. Each must exist.}
}
\value{
The diagnostic object with added columns ending in \code{_smooth} and a \code{smoothing} specification.
}
\details{
The maintained consolidated workflow uses raw exhaustive diagnostics.
This helper is not called by it, and consolidated summaries do not supply the
historical \code{meanH_B} column by default. A moving average across lambda
changes the displayed diagnostic and must not be described as the cell's
exhaustive agreement. For current filtering use \code{CECselectStableLambdas()}.
}
\seealso{\link{CECmigration}, \link{CECselectStableLambdas}.}
\examples{
old <- list(summary = data.frame(lambda = c(0.6, 1, 1.4),
  stab_ratio_mean = c(0.7, 0.95, 1), stab0_mean = c(0.8, 0.9, 1),
  stabB_mean = c(0.56, 0.855, 1), meanH_B = c(1.2, 1.1, 1)))
CECaddSmoothedDiagnostics(old)$summary
}
