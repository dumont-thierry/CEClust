\name{bands_partitions}
\alias{bands_partitions}
\title{List Distinct Retained Partitions and Their Diagnostic Witnesses}
\usage{bands_partitions(fit)}
\arguments{\item{fit}{A \code{cec_bands} result with its allocation archive.}}
\description{
Return one row for each distinct reference partition that has retained cells
in the final diagram. Exact partition identity, not group count alone, defines
duplicates.
}
\details{
If a partition has several retained regions or successful witnesses, the table
reports one recorded witness. Inspect \code{fit$tiles} and \code{fit$points} for
the full coverage. Partition names \code{P1}, \code{P2}, and so on are local to
this returned table and are not stable identifiers across different searches.
The \code{key} column is the archived candidate identifier.

The reported parameter pair, stability values and \code{witness_saturation}
belong to the diagnostic witness. Its saturation can be larger than the final display
threshold when diagnostic measurement permits low positive saturation. The
retained cells themselves always satisfy the display saturation threshold;
their range is reported separately in \code{min_retained_saturation} and
\code{max_retained_saturation}.
}
\value{A data frame containing the local name, candidate key, group count,
decreasing group sizes, witness ID, witness \code{C} and \code{lambda},
algorithmic and bootstrap stability, RSI, \code{witness_saturation},
\code{min_retained_saturation}, \code{max_retained_saturation}, and number of
retained cells. An empty data frame with these columns is returned if none
are retained.}
\seealso{\code{\link{bands_plot}}, \code{\link{bands_at}}}
