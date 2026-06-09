\name{CECclassif}
\alias{CECclassif}
\title{Fit one CEClust model}
\usage{
CECclassif(Z, lambda = 1, C = 1, Cquali = Inf, r0 = NULL,
  Nshots = 100, Nloop = 1000,
  familyType = "gaussAndDiscreteVector", sizeMaxOutlier = 0,
  autoRegroupOutliers = FALSE, displayRemainingTime = FALSE,
  focus = NULL, backend_data = NULL, n_cores = 1L,
  shot_strategy = c("current", "independent"))
}
\arguments{
\item{Z}{Data to cluster.}
\item{lambda}{Entropy regularisation parameter.}
\item{C}{Density bound for quantitative reference functions.}
\item{Cquali}{Optional frequency floor for qualitative variables.}
\item{r0}{Initial upper bound on the number of clusters.}
\item{Nshots}{Number of random starts.}
\item{Nloop}{Maximum number of descent iterations per start.}
\item{familyType}{Reference family. Common values are
\code{"gaussUniv"}, \code{"gaussVector"}, \code{"discreteVector"}, and
\code{"gaussAndDiscreteVector"}.}
\item{sizeMaxOutlier, autoRegroupOutliers, displayRemainingTime, focus,
backend_data, n_cores, shot_strategy}{Advanced algorithmic options.}
}
\value{
A list containing the selected partition, fitted parameters, criterion value,
and realised number of occupied clusters \code{REO}.
}
\description{
Fits a composite entropy clustering model for one pair \code{(C, lambda)} using
multiple random starts.
}
\examples{
set.seed(1)
Z <- matrix(rnorm(80), ncol = 2)
fit <- CECclassif(
  Z = Z,
  C = 2,
  lambda = 1,
  r0 = 4,
  Nshots = 2,
  Nloop = 5,
  familyType = "gaussVector"
)
fit$REO
}
