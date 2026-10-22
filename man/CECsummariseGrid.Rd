% Maintained by hand; see CONTRIBUTING.md.
\name{CECsummariseGrid}
\alias{CECsummariseGrid}
\title{Read the cell-by-cell grid results}
\description{
Return the stored summary without refitting, filtering again or choosing a unique partition.
}
\usage{
CECsummariseGrid(x, ...)
}
\arguments{
\item{x}{A grid returned by \code{CECfitBoundGrid()} or \code{CECfitPreset()}. For a consolidated one-C diagnostic, use its \code{grid} component.}
\item{...}{Reserved; currently unused.}
}
\value{
A data frame with one row per stored (C,lambda) pair. Key consolidated columns are listed in \link{CECresults}.
}
\details{
\code{keep} records the filtering decision made when the grid was built.
Changing thresholds in a plot or in \code{CECselectStableLambdas()} does not
alter this table. \code{sat} is a fraction in [0,1], not a percentage.
\code{partition_id} identifies identical labelings within the same grid object;
these identifiers need not match those of another fit.
}
\seealso{\link{CECresults}, \link{CECselectStableLambdas}, \link{CECextractPartition}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
g <- CECfitBoundGrid(Z = X, C_grid = c(0.2, 1),
  lambda_grid = c(0.6, 1, 1.5), optimization = "em_bank",
  penalty = "none", k0 = 3, B = 3, algo_seed = 14,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  n_cores = 1)
s <- CECsummariseGrid(g)
s[, c("C", "lambda", "REO", "sat", "measured", "keep")]
unique(s$partition_id[s$keep])
}
