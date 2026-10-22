% Maintained by hand; see CONTRIBUTING.md.
\name{CECresults}
\alias{CECresults}
\title{Read and preserve a consolidated grid object}
\description{
Object structure and recommended accessors for grids fitted with explicit optimization choices.
}
\section{Summary columns}{
\describe{
\item{\code{C, lambda}}{Cell coordinates.}
\item{\code{REO}}{Realized number of occupied groups.}
\item{\code{sat, n_sat_clusters, any_bound_saturation}}{Observation-weighted saturation, active-group count and activation flag.}
\item{\code{objective, entropy, penalty_value, Hphi}}{Hard penalized score, unpenalized H_lambda, relative penalty, and objective/lambda.}
\item{\code{stab_algo, stab_boot, rsi, measured}}{Exhaustive diagnostic values and availability flag.}
\item{\code{keep_stability, exclude_saturation, keep, fail_reason}}{Stored filtering decision and explanation.}
\item{\code{partition_id}}{Index of the representative in parts. Equal IDs within one object identify the same stored partition, up to its canonical labels.}
}
Use \code{CECsummariseGrid()} rather than reconstructing this table.
}
\section{Stored components}{
\describe{
\item{\code{Z}}{Original numeric fitting data, row order preserved.}
\item{\code{display_data}}{Display coordinates and optional external labels.}
\item{\code{banks}}{L+B independent execution records, original runs first. EM banks retain candidates and start records; fixed-grid CEM banks retain cell fits. Staged adaptive CEM uses \code{measurements} instead and leaves \code{banks} empty.}
\item{\code{parts}}{Pooled deduplicated projected partitions.}
\item{\code{winners}}{One row per summary cell, one column per independent execution; entries index parts. These generate diagnostics. In adaptive CEM an unmeasured row contains NA; no pilot result is inserted into a diagnostic column.}
\item{\code{selected}}{One pooled representative ID per summary cell. These generate the displayed diagram.}
\item{\code{bootstrap_indices}}{The B resample index vectors.}
\item{\code{runs_by_C}}{Display-compatible one-C results with best_parts and lambda_diag.}
\item{\code{meta}}{Optimizer, penalty, controls, seeds, thresholds, measurement policy, data fingerprint, package version and fitting time.}
}
Prefer \code{CECextractPartition()}, \code{CECextractBestPartitions()} and
the plotting functions for routine use. Internal records are useful for audits,
but should not be edited before rescoring.
}
\section{Adaptive results}{
Adaptive objects additionally inherit from \code{CEC_adaptive_grid}. The
\code{adaptive} component records the initial window, actual leaf geometry,
extension history, right-boundary evidence, refinement, harmonization controls,
optional quadrant audits and stopping limits. See \link{CECadaptiveDiagram}.
\code{meta$grid_strategy} distinguishes fixed and adaptive evaluations.
Adaptive CEM records \code{meta$protocol="staged_cem_v1"}, the pilot in
\code{adaptive$pilot} and \code{adaptive$pilot_fits}, and independent campaigns
in \code{measurements}, keyed by cell ID. Each campaign retains its runs,
bootstrap projection, best partition, diagnostics and saturation at the time
it was requested. \code{executions} lists phase, seed, cell, start use and
checkpoint reuse. \code{adaptive$final_refinement_pending} reports whether
measurement changes reactivate a frontier. Such a result is not marked as
completed adaptation.
}
\section{Saving and comparing}{
Use \code{saveRDS(g,"cec_none.rds")} and distinct names for each scoring
configuration; recover with \code{readRDS()}. Store the source data transformation,
package version and \code{sessionInfo()} alongside the object.
CECcompactGrid preserves the banks needed by CECrescoreGrid.
IDs are local to an object: compare labels using optimal agreement rather than
equating partition_id across separate runs. The internal fitting_elapsed field
includes fitting/checkpoint loading and, for staged CEM, harmonization; it
excludes agreement comparisons and does not include the whole scoring and
plotting pipeline. Staged CEM also records pilot_elapsed, measurement_elapsed,
diagnostics_elapsed and elapsed; these overlapping phase totals must not be
added together. Use
\code{system.time()} around the full call for an elapsed-time comparison.
}
\seealso{\link{CECsummariseGrid}, \link{CECextractPartition}, \link{CECrescoreGrid}, \link{CECcompactGrid}.}
