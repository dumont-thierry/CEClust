\name{bands_score}
\alias{bands_score}
\title{Evaluate the Profiled Gaussian Composite Entropy of a Hard Partition}
\usage{
bands_score(X, labels, C, lambda = 1, penalty = c("none", "bic"))
}
\arguments{
\item{X}{A finite numeric matrix, or a numeric vector for one dimension.}
\item{labels}{One group label per observation. Arbitrary renumbering does not
change the partition or the score.}
\item{C}{Positive density ceiling. \code{Inf} requests the unbounded regular
Gaussian profile.}
\item{lambda}{Finite nonnegative weight of the entropy of group proportions.}
\item{penalty}{\code{"none"} or the Gaussian order correction \code{"bic"}.}
}
\description{
Fit the Gaussian group parameters from the supplied hard allocations and
evaluate their profiled composite entropy. Smaller values are preferred.
}
\details{
The group covariance uses its maximum-likelihood divisor, equal to the group
size. For finite \code{C}, a regular covariance whose density peak exceeds the
ceiling is inflated proportionally until the ceiling is reached. Singular
groups contribute \eqn{-\log(C)} and are entirely saturated. No jitter or
density floor is added.

The score is the weighted Gaussian cross-entropy plus \code{lambda} times the
entropy of the group proportions. The optional order correction is
\eqn{(K-1)(d+1)(d+2)\log(n)/(4n)}, where \eqn{K} is the number of nonempty
groups. At an infinite ceiling a singular Gaussian group does not have a
finite unbounded profile; such a group must not be treated as a regular
unconstrained fit.
}
\value{A list with \code{J} (the criterion), \code{entropy}, \code{K},
\code{saturation}, and the additive \code{penalty}. Saturation is the total
mass of groups reaching the bound, not the proportion of groups.}
\examples{
x <- c(-2, -1, 0, 1, 2)
bands_score(x, rep(1, length(x)), C = 10)
bands_score(x, seq_along(x), C = 10, lambda = 1.3)
}
\seealso{\code{\link{cec_bands}}}
