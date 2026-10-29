\name{bands_control}
\alias{bands_control}
\title{Control Candidate Search, Diagram Refinement and Diagnostics}
\usage{
bands_control(orders = 1:20, starts = 5L, extra_starts = 5L,
  max_order = 320L, max_iter = 1000L, tolerance = 1e-8,
  max_anchors = 40L, C_points = 41L, max_C_points = 257L,
  refine_rounds = 2L, projection_rounds = 2L,
  max_projection_queries = 2000L, cache_mb = 64,
  diag_mode = c("components", "all"), algo_min = .85,
  rsi_min = .95, saturation_max = 0,
  measurement_saturation = .10, branches = c("free", "bounded"))
}
\arguments{
\item{orders}{Positive integer initial group counts for the pilot search.}
\item{starts}{Number of starts per pilot order.}
\item{extra_starts}{Additional random starts per anchor on each analysed sample.
Each start is used by both selected descent branches. Zero disables these
additional starts.}
\item{max_order}{Maximum requested initial group count, including the doubled
anchor orders. Actual starts also respect the available observations. This is
a search budget, not a bound on the analytically represented singleton endpoint.}
\item{max_iter}{Maximum number of CEM iterations per descent.}
\item{tolerance}{Positive relative tolerance used when accepting a decrease of
the objective.}
\item{max_anchors}{Maximum number of pilot anchors used by each repeated search.
When this budget is binding, the limitation is reported.}
\item{C_points}{Initial number of logarithmic density-ceiling values.}
\item{max_C_points}{Maximum number of density-ceiling values after inserting
critical scales and refining neighbouring columns.}
\item{refine_rounds}{Maximum number of rounds refining changes of partition or
saturation between columns. Zero disables these rounds.}
\item{projection_rounds}{Maximum number of bootstrap-projection enrichment
rounds of the pooled reference diagram. Zero disables these rounds.}
\item{max_projection_queries}{Budget for bootstrap-projection queries used to
enrich the reference collection. One query is one bootstrap bank at one
parameter pair. All bootstrap banks at a point are processed, or none if the
remaining budget is insufficient. A binding budget is a search limitation.}
\item{cache_mb}{Positive size budget in MiB for a transient cache. Detailed
allocation archives remain on disk rather than in an unlimited cache.}
\item{diag_mode}{\code{"components"} tests eligible cells from the barycentre
outward in each connected region and stops at the first successful witness.
\code{"all"} tests every eligible final cell.}
\item{algo_min, rsi_min}{Strict lower thresholds for algorithmic stability and RSI.}
\item{saturation_max}{Largest saturation allowed in retained cells; the default
is exactly zero.}
\item{measurement_saturation}{Largest saturation allowed when selecting a cell
for diagnostic evaluation. Must be at least \code{saturation_max}.}
\item{branches}{Descents included in candidate generation. The default uses
both \code{"free"} and \code{"bounded"}; a single branch can be selected for
an explicit comparison of search procedures.}
}
\description{
Create a checked control list for \code{\link{cec_bands}}. Search budgets,
measurement coverage and diagnostic thresholds are separate settings.
}
\details{
The free branch uses no finite density ceiling. Singular initial groups are
invalid for that branch; a singular proposal is not accepted as a finite
unconstrained optimum. The bounded branch evaluates the finite density ceiling
at the anchor. Neither branch introduces a hidden covariance ridge.

The diagnostic procedure has no champion-agreement rejection and no conditional
post-hoc replacement of RSI. In component mode, a successful witness supports
a regional display decision, not a claim that all cells were individually
measured. Saturation is still filtered at each displayed cell.
}
\value{A list containing the checked settings.}
\examples{
bands_control(orders = 1:5, starts = 2, extra_starts = 2,
              max_order = 10, max_anchors = 8)
bands_control(diag_mode = "all")
}
\seealso{\code{\link{cec_bands}}}
