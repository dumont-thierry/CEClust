\name{cec_bands}
\alias{cec_bands}
\title{Build a Gaussian CEC Diagram with Symmetric Bootstrap-Anchored Searches}
\usage{
cec_bands(X, L = 20L, B = 20L, workers = 1L, seed = 1L,
  bootstrap_indices = NULL, pilot_indices = NULL,
  control = bands_control(), penalty = c("none", "bic"),
  directory = tempfile("bands-"), resume = TRUE)
}
\arguments{
\item{X}{A finite numeric matrix, or a numeric vector for one dimension.
Coordinates and row order are used exactly as supplied.}
\item{L}{Number of independent original-data banks, at least two.}
\item{B}{Number of independent outer-bootstrap banks, at least one.}
\item{workers}{Maximum number of parallel bank workers. Diagram construction
and diagnostic assembly are separate phases.}
\item{seed}{Nonnegative integer master seed. Separate deterministic seeds are
derived for the pilot, each bank, the internal bootstraps and additional starts.
The caller's R random-number state is restored on exit.}
\item{bootstrap_indices}{Optional list of exactly \code{B} integer index vectors,
each of length \code{nrow(X)}, with replacement permitted. These define the
outer bootstrap samples, not the internal initialization bootstraps.}
\item{pilot_indices}{Optional index vector of length \code{nrow(X)} defining the
pilot bootstrap. Otherwise it is drawn using its own derived seed.}
\item{control}{A control list returned by \code{\link{bands_control}}.}
\item{penalty}{\code{"none"} or the Gaussian BIC-type order correction.}
\item{directory}{Analysis directory for candidate allocations, checkpoints,
configuration, tables and the progress page. Keep this directory with the
returned object.}
\item{resume}{If true, reuse compatible completed checkpoints in the analysis
directory. Incompatible inputs, controls, seeds or runtime are rejected. If
false, an existing analysis manifest is rejected; it is not overwritten.}
}
\description{
Generate a pilot on an initial bootstrap, select anchors in its winning bands,
and run the same randomized estimator on original and outer-bootstrap samples.
Candidates are compared through a pooled reference diagram while each
repetition retains its own bank for reproducibility diagnostics.
}
\details{
For every analysed sample, an internal bootstrap perturbs each pilot allocation
rule. Free and bounded CEM descents are run from the same initial allocations
on that internal bootstrap, followed by projection and descent on the analysed
sample. Additional starts at twice the pilot group count receive the same two
branches. Finite order and sample-capacity restrictions are recorded. Exact
partition duplicates are removed independently of arbitrary group numbering.

The final reference collection combines original-bank candidates with outer-
bootstrap winners projected onto the original observations at recorded anchor
and diagram queries. Projected hard partitions have their parameters refitted
on the original sample. This enrichment has explicit round and query budgets;
it does not enumerate every possible projected allocation on a continuous
parameter domain. The final adaptive mesh is constructed from the reference
collection, rather than frozen to the pilot mesh.

At each measured point, algorithmic stability uses all \eqn{L(L-1)/2}
original-original agreements and bootstrap stability uses all \eqn{LB}
original-bootstrap agreements. RSI is their ratio and can exceed one.
Agreement maximizes matched observations over one-to-one group-label matchings.
There is no champion-agreement gate or conditional post-hoc RSI replacement.

The common pilot is held fixed. Independent-run diagnostics are conditional on
that pilot and do not measure the variability of relearning it. The pooled
champion is selected using the original and bootstrap banks, so these same
diagnostics are not an independent validation of its selection.

Candidate rules and scores are held compactly; hard allocations and audit
checkpoints are archived on disk. A transient bounded cache is excluded from
the returned scientific object. The progress page is \code{SUIVI.html} in
\code{directory}. The bands engine is included in CEClust.
}
\value{
An object of class \code{cec_bands}. Its components include \code{X},
\code{control}, \code{pilot}, \code{anchors}, independent \code{runs}, the pooled
\code{reference}, final \code{tiles}, measured \code{points}, \code{store},
\code{timings}, and explicit search \code{limits}. Each bank records its outer
and internal indices, seed, candidate catalogue and descent log. Detailed
point checkpoints contain the independent winners and all pairwise agreements.
See \code{\link{bands_partitions}} for a compact table of retained partitions.
}
\examples{
\dontrun{
set.seed(12)
x <- c(rnorm(60, -2, .6), rnorm(60, 2, .6))
control <- bands_control(orders = 1:4, starts = 2, extra_starts = 2,
                         max_order = 8, max_anchors = 6,
                         C_points = 11, refine_rounds = 1,
                         projection_rounds = 1,
                         max_projection_queries = 100)
fit <- cec_bands(x, L = 4, B = 4, control = control,
                 directory = "bands-example")
bands_plot(fit)
bands_partitions(fit)
bands_validate(fit)
}
}
\seealso{\code{\link{bands_control}}, \code{\link{bands_at}},
\code{\link{bands_plot}}, \code{\link{bands_validate}}}
