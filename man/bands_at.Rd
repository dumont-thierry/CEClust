\name{bands_at}
\alias{bands_at}
\title{Evaluate Independent-Bank Diagnostics at a Specified Parameter Pair}
\usage{bands_at(fit, C, lambda)}
\arguments{
\item{fit}{A \code{cec_bands} result whose analysis archive remains available.}
\item{C, lambda}{Positive finite scalar values of the density ceiling and
entropy weight.}
}
\description{
Select the winner of each independent bank at the supplied parameter pair,
project bootstrap rules onto the original observations, and calculate all
original-original and original-bootstrap agreements.
}
\details{
No CEM descent is launched and the frozen reference collection is not enriched.
Queries outside the displayed mesh still use the same finite candidate banks;
they do not prove coverage of undiscovered partitions there. Newly computed
allocations may be written to the analysis archive for traceability.

The \code{passed} field tests the stability and RSI thresholds only. The
reference saturation is reported separately; retention additionally requires
the requested saturation condition. A bootstrap rule with unsupported original
observations makes the diagnostic invalid rather than silently reallocating
those observations using regularization.
}
\value{A list containing validity, threshold passage, reference identity and
score, reference saturation, independent winners, and the exhaustive diagnostics.
Invalid outputs include a reason and the failed run. No conditional RSI
replacement is performed.}
\seealso{\code{\link{cec_bands}}, \code{\link{bands_partitions}}}
