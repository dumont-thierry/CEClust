\name{bands_cleanup}
\alias{bands_cleanup}
\title{Remove retained partitions challenged by group subdivisions}
\description{
Apply a separate split-and-CEM removal filter to the partitions retained in a
CEC diagram. The original fit, independent banks and stability diagnostics are
unchanged. No diagram reconstruction or replacement partitions are introduced.
}
\usage{
bands_cleanup(fit, directory, workers = 1L, max_rounds = 4L,
  trials = 4L, seed = 1L, score_tolerance = 1e-8,
  entropy_tolerance = 1e-8, resume = TRUE,
  reuse_directory = NULL)
}
\arguments{
\item{fit}{A result from \code{cec_bands}, with its archive available.}
\item{directory}{A separate directory for the filter's results and checkpoints.}
\item{workers}{Number of parallel CEM workers; at most \code{trials} jobs run concurrently.}
\item{max_rounds}{Maximum recursive split rounds per partition.}
\item{trials}{At least two trials per round. The last two split all splittable
groups. Preceding trials split one random group, then random halves. Odd trials
use a median cut along the leading within-group principal direction; even trials
use two-means with a median-cut fallback.}
\item{seed}{Nonnegative integer controlling group and split choices.}
\item{score_tolerance}{Positive relative numerical margin for score improvement.}
\item{entropy_tolerance}{Nonnegative minimum increase in allocation entropy.}
\item{resume}{Reuse checkpoints when the fit, runtime, search settings and seeds
match. Otherwise an existing cleanup archive is rejected.}
\item{reuse_directory}{Optional earlier cleanup archive directories. Source archive identity,
source runtime, numerical settings, penalty and seed must match. A descent is
reused only when its exact split initialisation and parameter pair also match.}
}
\details{
Each retained partition is tested at its original representative pair
\eqn{(C,\lambda)}. Selected groups are split in two, then native CEM is run with
the fit's criterion, convergence settings and order penalty. Recursive split
seeds are retained even when CEM merges the children, so later rounds can test
higher orders. The original bank's \code{max_order} does not cap these seeds.

A partition is removed as soon as a completed trial has a criterion lower by
more than \code{score_tolerance * max(1, abs(original_score))} and either its
allocation entropy increases by more than \code{entropy_tolerance}, or its
saturation exceeds \code{fit$control$saturation_max}. With the usual
zero-saturation retention limit, the latter means positive saturation.
The score improvement is mandatory: saturation alone with a worse score never
suffices. A score improvement without higher entropy or saturation does not
suffice either. Each condition is stored separately, together with score gain,
competitor order and saturation.

This is a filter at the existing representative point. It does not establish
that a removed partition loses throughout its old parameter region. Surviving
the finite split budget is not an optimality certificate. The function tests
all retained partitions, not every losing candidate in the bank.

The original result is not edited, no replacement partition is added, and no
new stability or RSI is calculated. The diagnostic columns in the returned
table are the original diagnostics. The filter is an explicit post-fit call,
not automatically run by \code{cec_bands}. Its elapsed time is additional to
the original fit. Reuse or restart does not constitute a fresh timing experiment.
}
\value{
A \code{bands_cleanup} list containing \code{partitions} (old partitions and
removal decisions), \code{remaining} (only those not removed), \code{descents},
per-partition \code{targets}, \code{limits}, archive paths and elapsed time.
\code{cleanup_result.rds}, \code{cleanup_partitions.csv}, and
\code{cleanup_descents.csv} are saved in \code{directory}. Every trial has a
checkpoint. Keep the original fit archive together with the cleanup archive.
}
\examples{
\dontrun{
audit <- bands_cleanup(fit, "my-analysis-cleanup", workers = 4,
                      max_rounds = 4, trials = 4, seed = 123)
audit$partitions[, c("name", "K", "removed", "reason")]
audit$remaining
}
}
\seealso{\code{\link{cec_bands}}, \code{\link{bands_partitions}}}
