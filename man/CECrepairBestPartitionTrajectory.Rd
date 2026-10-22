% Maintained by hand; see CONTRIBUTING.md.
\name{CECrepairBestPartitionTrajectory}
\alias{CECrepairBestPartitionTrajectory}
\title{Repair a historical representative trajectory}
\description{
Compatibility operation that propagates better historical unpenalized candidates along a lambda path. It is not an additional stage of the EM-bank workflow.
}
\usage{
CECrepairBestPartitionTrajectory(best_partitions_obj, Z, tol = 1e-10, max_iter = 100)
}
\arguments{
\item{best_partitions_obj}{Representative path from \code{CECextractBestPartitions()}.}
\item{Z}{Original fitting observations, in their original order.}
\item{tol}{Score comparison tolerance.}
\item{max_iter}{Maximum number of historical propagation passes.}
}
\value{
The representative object with repair metadata. An EM-bank object is returned with zero repairs and a reason indicating minimization over the common bank.
}
\details{
All EM-bank cells already search the same candidate pool; this repair
is unnecessary there. Penalized CEM paths are rejected because the historical
repair does not evaluate their penalized objective. On an unpenalized historical
path, repair changes representatives only: it does not retrospectively recompute
independent execution diagnostics.
}
\seealso{\link{CECrescoreGrid}, \link{CECextractBestPartitions}, \link{CECmigration}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
d <- CECfitLambdaGrid(X, C = 1, lambda_grid = c(0.6, 1, 1.5),
  optimization = "em_bank", penalty = "none", k0 = 3, B = 3,
  seed = 14, optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  verbose = FALSE)
best <- CECextractBestPartitions(d)
repaired <- CECrepairBestPartitionTrajectory(best, X)
repaired$repair_reason
}
