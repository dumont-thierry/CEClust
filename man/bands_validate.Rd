\name{bands_validate}
\alias{bands_validate}
\title{Check Archived Allocations and Recalculate the Stored Agreements}
\usage{bands_validate(fit)}
\arguments{\item{fit}{A \code{cec_bands} result with its complete analysis archive.}}
\description{
Check bank roles and seeds, internal bootstrap indices, allocation checksums,
retained-cell saturation, and agreement counts. Recalculate stored diagnostics
from archived allocations with caching disabled.
}
\details{
The validation performs no new CEM fitting. It reads archived allocations and
point checkpoints, so its cost increases with the number of stored candidates
and measured points. Run it outside a timed comparison with other methods.

This is an integrity and arithmetic check. A passing result does not prove
global optimality, exhaustive candidate coverage, or scientific relevance of
each retained partition. Missing or incompatible archives cause an error.
}
\value{On success, a list with \code{passed = TRUE}, the numbers of banks,
measured points and retained partitions, and the recorded search limits.
An inconsistency raises an error.}
\seealso{\code{\link{cec_bands}}, \code{\link{bands_at}}}
