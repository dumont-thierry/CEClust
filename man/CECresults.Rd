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
\item{\code{banks}}{L+B independent execution records, original runs first. EM banks retain candidates and start records; CEM banks retain cell fits.}
\item{\code{parts}}{Pooled deduplicated projected partitions.}
\item{\code{winners}}{One row per summary cell, one column per independent execution; entries index parts. These generate diagnostics.}
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
}
\section{Saving and comparing}{
Use \code{saveRDS(g,"cec_none.rds")} and distinct names for each scoring
configuration; recover with \code{readRDS()}. Store the source data transformation,
package version and \code{sessionInfo()} alongside the object.
CECcompactGrid preserves the banks needed by CECrescoreGrid.
IDs are local to an object: compare labels using optimal agreement rather than
equating partition_id across separate runs. The internal fitting_elapsed field
does not include the entire scoring, diagnostic and plotting pipeline; use
\code{system.time()} around the full call for an elapsed-time comparison.
}
\seealso{\link{CECsummariseGrid}, \link{CECextractPartition}, \link{CECrescoreGrid}, \link{CECcompactGrid}.}
