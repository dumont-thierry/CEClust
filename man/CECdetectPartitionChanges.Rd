% Maintained by hand; see CONTRIBUTING.md.
\name{CECdetectPartitionChanges}
\alias{CECdetectPartitionChanges}
\title{Describe changes along a representative path}
\description{
Mark changes in group count, partition agreement, or both along a lambda path. This describes the fitted path and does not replace the stability filter.
}
\usage{
CECdetectPartitionChanges(best_partitions_obj, partition_metric = c("match", "ARI"), sim_threshold = 0.8, 
    criterion = c("REO_or_threshold", "REO_only", "threshold_only", "REO_and_threshold"), 
    include_first = TRUE, method = c("consecutive", "regime_reference"), lambda_subset = NULL)
}
\arguments{
\item{best_partitions_obj}{Output of \code{CECextractBestPartitions()}.}
\item{partition_metric}{\code{"match"} or \code{"ARI"}; see \link{CECstability}.}
\item{sim_threshold}{An agreement strictly below this threshold counts as a partition change.}
\item{criterion}{\code{"REO_only"}, \code{"threshold_only"}, their logical union \code{"REO_or_threshold"}, or their conjunction \code{"REO_and_threshold"}.}
\item{include_first}{Mark the first path entry as a change/reference.}
\item{method}{\code{"consecutive"} compares adjacent entries. \code{"regime_reference"} compares with the start of the current regime and updates that reference when a change is detected.}
\item{lambda_subset}{Optional exact subset of stored lambda values to inspect.}
}
\value{
A list with a per-lambda \code{summary}, \code{change_lambdas}, \code{change_indices} and the chosen rules.
}
\details{
Use a path with increasing lambda. Equal REO does not imply equal
partitions; an agreement-based criterion can detect a regrouping without an
order change. Conversely, a display regime based on an agreement threshold
is not a collection of exactly identical partitions.
}
\seealso{\link{CECextractBestPartitions}, \link{CECselectStableLambdas}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
d <- CECfitLambdaGrid(X, C = 1, lambda_grid = c(0.6, 1, 1.5),
  optimization = "em_bank", penalty = "none", k0 = 3, B = 3,
  seed = 14, optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  verbose = FALSE)
changes <- CECdetectPartitionChanges(CECextractBestPartitions(d),
  criterion = "REO_or_threshold", sim_threshold = 0.95)
changes$change_lambdas
}
