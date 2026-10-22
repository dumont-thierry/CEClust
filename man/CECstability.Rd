% Maintained by hand; see CONTRIBUTING.md.
\name{CECstability}
\alias{CECstability}
\title{Compute exhaustive agreement and relative stability}
\description{
Calculate algorithmic stability, bootstrap stability and their ratio from partitions of the same reference observations, independently of the clustering method that produced them.
}
\usage{
CECstability(original, bootstrap, partition_metric = c("match", "ARI"))
}
\arguments{
\item{original}{List of L label vectors from independent executions on the original sample.}
\item{bootstrap}{List of B label vectors obtained by fitting bootstrap samples and projecting their outputs onto the same original observations. May be an empty list.}
\item{partition_metric}{\code{"match"}: maximum fraction of agreeing labels after optimal one-to-one relabeling, padding unequal group counts with empty labels. \code{"ARI"}: adjusted Rand index; it can be negative.}
}
\value{
A list with \code{Stab_algo}, \code{Stab_boot}, \code{RSI},
\code{original_pairs}, \code{cross_pairs}, \code{algorithm_agreements},
\code{bootstrap_agreements} and \code{definition}.
}
\details{
Every vector must be nonempty, complete and label the same observations
in the same order. Group numbers need not agree between runs. This function
does not fit, resample or project data.

Algorithmic stability averages the L(L-1)/2 distinct original-original pairs.
Bootstrap stability averages the L B original-bootstrap pairs. There are no
self-comparisons and no reference selected by score. RSI is Stab_boot/Stab_algo
when the denominator is finite and positive. Missing pairs give NA; in particular
L=1 cannot establish algorithmic stability, and B=0 cannot establish bootstrap
stability. The ratio is not clipped to one. See \link{CECdiagnostics} for interpretation.
}
\seealso{\link{CECdiagnostics}, \link{CECfitBoundGrid}.}
\examples{
original <- list(c(1, 1, 2, 2), c(4, 4, 3, 3), c(1, 2, 2, 2))
bootstrap <- list(c(1, 1, 2, 2), c(1, 1, 1, 2))
CECstability(original, bootstrap)
CECstability(original[1], bootstrap)$RSI  # NA: only one original run
}
