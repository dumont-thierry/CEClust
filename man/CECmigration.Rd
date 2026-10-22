% Maintained by hand; see CONTRIBUTING.md.
\name{CECmigration}
\alias{CECmigration}
\title{Use explicit solvers while retaining historical workflows}
\description{
Compatibility differences between historical linked-path calculations and the consolidated numeric engines in CEClust 0.3.0.9003.
}
\section{Choosing the route}{
For new numeric analyses, explicitly pass optimization="cem" or
"em_bank" and penalty="none" or "bic" to a grid function. Omitting all new
choices and using an empty optimizer_control retains the historical linked-path
engine. Explicit penalty or RSI policy also selects the consolidated route.
For \code{CECfitBoundGrid}, an explicit \code{grid_strategy} or a nonempty
\code{adaptive_control} also selects the consolidated engine.
The single-fit CEM function has an additional compatibility rule: unpenalized
CEM with no control except an optional seed keeps historical starts; a further
numeric control selects the consolidated CEM search.
}
\section{What changes}{
The maintained grid route uses independent repetitions and exhaustive
original-original/original-bootstrap comparisons. Its defaults are 20/20
repetitions and filtering thresholds 0.85/0.95/0. The historical full-grid route
has different budgets, linked paths and optional repair stages; the one-C
interface has its own explicit defaults. The historical diagnostic default
has also changed to exhaustive. Old saved values are not silently recalculated.
Legacy pairwise and reference_bestH definitions remain explicitly named.

The new EM bank is evaluated on the whole supplied grid without a final
harmonization. Reuse it with CECrescoreGrid. Historical warm-start, final-descent
and smoothing instructions should not be copied into that workflow. Unsupported
active options are rejected; some historical scheduling and display switches
are accepted without affecting the numeric engine.
}
\section{Compatibility functions}{
CECaddSmoothedDiagnostics and CECidentifyStableLambdas retain old
descriptive rules; current filtering uses CECselectStableLambdas.
CECextractBestPartitions requires source="all", criterion="projected_H"
for solver diagnostics. CECrepairBestPartitionTrajectory does nothing to an
EM-bank path and rejects penalized CEM paths.
The plot-path function still has historical default thresholds for its optional
one-dimensional allocation view: pass an explicit lambda_stable_obj obtained
from CECselectStableLambdas to use the grid's current filter.
}
\section{Adaptive CEM correction in 0.3.0.9003}{
Versions 0.3.0.9001 and 0.3.0.9002 fitted L+B repetitions at every evaluated
adaptive CEM cell, including pilot parents and cells outside the measurement
policy. Version 0.3.0.9003 replaces that orchestration with a separate pilot
and deferred selective campaigns. It restores parent starts in the pilot,
post-measurement improvements, local harmonization and coverage checks.
The default adaptive CEM policy is now saturation_band. Fixed-grid and EM-bank
protocols retain their defaults. Old adaptive CEM checkpoints are not accepted
as staged campaigns and old saved results are not silently recomputed.
See \link{CECadaptiveDiagram} for the stages and stored records.
}
\section{Scope of this version}{
New solvers are numeric Gaussian implementations. The package's historical
unpenalized CEM still supports categorical and mixed observations. Choose \code{grid_strategy="adaptive"} to use the integrated extension, right-row stopping and frontier refinement rules. Fixed-grid calls preserve their previous behaviour. See \link{CECadaptiveDiagram}.
Only the code in this package is required for its new engines; no local research
script or external precompiled EM DLL must be sourced.
}
\seealso{\link{CECfitBoundGrid}, \link{CECfitLambdaGrid}, \link{CECoptimizers}.}
