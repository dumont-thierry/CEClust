\name{CECfitLambdaGrid}
\alias{CECfitLambdaGrid}
\title{Fit a CEClust lambda grid}
\usage{
CECfitLambdaGrid(Z, lambda_grid = seq(0.1, 2, by = 0.1),
  k0 = 10, B = 10, C = 1, Cquali = Inf, r0 = NULL,
  Nshots_fresh = 5, Nshots_warm = 5, Nshots_init = NULL,
  Nloop = 300, familyType = "gaussAndDiscreteVector",
  sizeMaxOutlier = 0, autoRegroupOutliers = FALSE, focus = NULL,
  replace = TRUE, seed = NULL, silent = TRUE, verbose = TRUE,
  n_cores = 1L, checkpoint_dir = NULL, auto_checkpoint = NULL,
  resume = TRUE, batch_size = NULL, show_progress = NULL,
  rerun_failed_serial = TRUE, checkpoint_compress = FALSE,
  partition_metric = c("match", "ARI"),
  stability_approach = c("pairwise", "reference_bestH"))
}
\arguments{
\item{Z}{Data to cluster.}
\item{lambda_grid}{Values of \code{lambda} to evaluate.}
\item{k0}{Number of linked trajectories on the original data.}
\item{B}{Number of bootstrap trajectories.}
\item{C}{Density bound.}
\item{Cquali}{Optional qualitative frequency floor.}
\item{r0}{Initial upper bound on the number of clusters.}
\item{Nshots_fresh, Nshots_warm, Nshots_init}{Numbers of random or warm starts
used along the lambda path.}
\item{Nloop}{Maximum number of descent iterations per local fit.}
\item{familyType}{Reference family.}
\item{sizeMaxOutlier, autoRegroupOutliers, focus, replace, seed, silent,
verbose, n_cores, checkpoint_dir, auto_checkpoint, resume, batch_size,
show_progress, rerun_failed_serial, checkpoint_compress, partition_metric,
stability_approach}{Advanced workflow options.}
}
\value{
A diagnostic object containing fitted lambda-path partitions, stability
diagnostics, RSI, saturation diagnostics, and summaries.
}
\description{
Fits CEClust along a grid of \code{lambda} values for one fixed value of
\code{C}. For complete \code{C x lambda} studies, use
\code{\link[=CECfitBoundGrid]{CECfitBoundGrid()}}.
}
\examples{
set.seed(1)
Z <- matrix(rnorm(80), ncol = 2)
lambda_grid <- CECfitLambdaGrid(
  Z = Z,
  lambda_grid = c(0.8, 1.2),
  k0 = 1,
  B = 1,
  C = 2,
  r0 = 4,
  Nshots_fresh = 1,
  Nshots_warm = 1,
  Nloop = 5,
  familyType = "gaussVector",
  n_cores = 1,
  show_progress = FALSE
)
head(lambda_grid$summary)
}
