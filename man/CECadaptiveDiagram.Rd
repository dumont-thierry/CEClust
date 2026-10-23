\name{CECadaptiveDiagram}
\alias{CECadaptiveDiagram}
\title{Construct and read an adaptive CEC diagram}
\description{
Extend an initial window in four directions and refine the cells on both sides
of relevant boundaries. This explores a collection of partitions; it does not
select a single optimal pair of tuning parameters.
}
\section{Entry point}{
Use \code{CECfitBoundGrid(grid_strategy="adaptive", ...)}. Both
\code{optimization="cem"} and \code{"em_bank"} support \code{penalty="none"}
and \code{"bic"}. These adaptive engines require numeric Gaussian data.
No external research script is loaded.
\code{CECrescoreGrid()} preserves the adaptive strategy by default and can
reconstruct the mesh with a different penalty from the same EM banks.
}
\section{Initial window}{
Supply at least two positive C values equally spaced on the log scale and two
positive lambda values equally spaced on the linear scale. The lower lambda
cell edge must remain positive. The supplied coordinates are initial cell
centres, not limits on the final window.
If omitted, five C values range from one quarter to four times the fitted
one-group Gaussian peak, and five lambda values range from 0.5 to 2.
Singular data require an explicit initial C grid.
For EM, specify \code{optimizer_control$Cmax} before fitting. Its adaptive
default is 64 times the largest initial C. The bank is generated once under
this ceiling. Evaluated cell centres cannot exceed it; reaching it is a
reported limit, not evidence that no structure exists further right.
}
\section{Extension rules}{
Left extension stops when every boundary cell is trivial or exceeds the final
saturation threshold. Lower extension stops when every boundary cell exceeds
that threshold; upper extension stops when every boundary cell is trivial.
For the right boundary, each lambda row must be trivial or exceed the threshold,
or have zero saturation and unchanged REO over three successive evaluated
columns (two extensions). Rows must have exactly matching lambda intervals.
Refining the comparison columns or changing those intervals resets their history.
The REO rule is a practical stopping rule, not a guarantee that partitions agree
or that all structure beyond the window has been found.
}
\section{Refinement and stopping}{
Cells sharing a positive-length edge are neighbours. A frontier separates an
eligible and an ineligible saturation regime, or two different partitions
within the eligible regime. Both neighbouring cells are subdivided into four.
Different partitions with the same REO are distinguished by default.
Refinement becomes locally inactive after successive rounds introduce no new
local winner and move the frontier by less than the tolerance. A changed
partition or saturation regime wakes that neighbourhood again.
The displacement measure is numerical, not a bound on an unknown true frontier.
}
\section{Adaptive controls}{
Pass a named \code{adaptive_control} list. Its fields are:
\describe{
\item{\code{max_cells=12000}}{Maximum number of active rectangular cells.}
\item{\code{max_extensions=16}}{Maximum extensions per side; zero disables extension.}
\item{\code{max_rounds=12, max_depth=8}}{Refinement budgets. The initial mesh is round one.}
\item{\code{patience=2, position_tolerance=0.125}}{Local stopping controls; displacement is measured in initial-grid units.}
\item{\code{distinguish_partitions=TRUE}}{Refine partition changes even when REO is unchanged.}
\item{\code{frontier_saturation}}{Saturation threshold for relevant frontiers; defaults to sat_threshold and cannot be smaller.}
\item{\code{harmonize=TRUE, harmonization_min=0.95}}{For CEM, propagate improving neighbour partitions and run a descent from them. The initial pass visits the whole mesh. Later passes start from the newly evaluated cells and their neighbours; an improvement reactivates affected neighbours until propagation stops. Keep requires agreement at least harmonization_min between the best independent campaign partition and the final harmonized representative. It is missing for cells without a campaign. Independent repetitions remain unchanged. This additional agreement filter applies to the cell strategy unless the conditional reference branch is used. Regional acceptance uses the two diagnostic thresholds only. EM uses its common bank and has no additional harmonization.}
\item{\code{diagnostic_strategy="components"}}{Freeze the final pilot and traverse connected components of exactly the same partition. Use "cells" for the previous cell-by-cell strategy.}
\item{\code{reference_bootstrap=TRUE}}{Apply the conditional reference-bootstrap rule to every tested point before accepting it; see \link{CECdiagnostics}. This adaptive setting overrides the top-level value when supplied.}
\item{\code{trivial_quadrants=FALSE}}{Optional heuristic exemption from diagnostic measurement in validated upper-left trivial quadrants.}
\item{\code{audit_fraction=0.10, anchor_stability=0.95}}{Quadrant controls described below.}
}
Inspect \code{g$adaptive$history}, \code{margins}, \code{restrictions},
\code{stop_reason} and \code{completed}. A budget stop or a Cmax stop is never
reported as completed adaptation. Even completed adaptation only means that
the specified numerical stopping rules were met within the explored window.
}
\section{Connected-component traversal}{
The default \code{diagnostic_strategy="components"} freezes the harmonized pilot
before diagnostic campaigns. Edge-sharing cells belong to one component only
when their label vectors describe exactly the same partition, up to relabeling.
Disconnected regions remain separate, even when they share a partition.
Saturation does not split components.

The barycentre is computed over the whole component in (log C, lambda), weighting
each cell centre by its rectangular area in these coordinates. This prevents
unequal refinement from moving the geometric centre. Candidate cells are ordered
by increasing Euclidean distance to this barycentre, without rescaling the axes;
cell ID breaks ties. The measurement policy then restricts the candidates.
Thus the default CEM saturation band still excludes test points above 0.10.

Each candidate receives the ordinary exhaustive diagnostics and, when all three
conditions hold, the conditional reference-bootstrap recalculation. This is done
before comparing Stab_algo and RSI to their strict thresholds. A component stops
at its first passing point, or after exhausting its eligible candidates.
No extra A_harm filter is applied to regional acceptance; its direct value is
retained for interpretation. Campaign winners never alter the frozen pilot.

For display and selection, a passing witness supplies its diagnostic values to
the whole component. The final saturation mask is applied separately. This is
a regional approximation: it does not establish stability at each untested point.
\code{measured} identifies actual campaigns, \code{diagnostic_origin} identifies
regional values, and \code{diagnostic_witness} gives their measured source.
\code{stab_algo_direct}, \code{stab_boot_direct}, and \code{rsi_direct} preserve
pointwise values, including failed attempts. \code{adaptive$components} records
membership, full barycentres, candidate paths, attempts and stopping decisions.
\code{winners} retains only the independent outputs at measured CEM cells;
untested CEM rows remain missing. EM winners remain available from its banks.

This strategy and optional trivial-quadrant screening are alternatives and
cannot be combined. Tightening plotting thresholds cannot reconstruct skipped
campaigns; rerun diagnostics if complete coverage at new thresholds is required.
}
\section{Diagnostics and optional quadrants}{
With \code{diagnostic_strategy="cells"}, adaptive CEM has two stages. The pilot runs one complete multistart search
on the original observations at each new cell, followed by neighbour
harmonization. A child cell additionally starts from its parent's partition.
No original-data repetition campaign or bootstrap fit is run during this
calibration and refinement stage.

After the pilot stops, independent L original and B bootstrap executions are
run only at active cells satisfying the measurement policy. The default for
adaptive CEM is \code{"saturation_band"}, with a maximum of 0.10; explicit
\code{"all"} requests campaigns even in saturated cells. Campaigns have
separate seeds and receive no pilot, parent or neighbour starts. Bootstrap
models are projected onto the original observations with frozen parameters.
The two exhaustive agreement averages use these frozen independent outputs.

A better campaign partition can improve the representative diagram. Such
improvements are harmonized locally; borders and measurement coverage are
checked again. Newly eligible cells receive campaigns, while already measured
cells are never refitted or substituted by a common winner. Post-measurement
changes can reactivate a refinement boundary: inspect
\code{final_refinement_status} and \code{final_refinement_pending}.
A pending boundary or a budget stop prevents \code{completed=TRUE}.

EM uses a different arrangement: L+B independent banks are fitted once and
shared across cells. Their winners are recomputed as the grid changes; the
measurement policy gates agreement calculations, not bank generation. No
extra CEM harmonization is applied. Identical winner patterns and partition
comparisons are cached exactly.

Optional trivial quadrants use a measured, trivial anchor passing RSI and
algorithmic-stability thresholds. Every explored cell with C no greater and
lambda no smaller must be trivial. Controls use one deterministic, persistent
budget for the union of candidate quadrants; overlapping quadrants do not
duplicate that budget. A failing control invalidates its quadrant. An exempted
cell retains missing diagnostics and is marked \code{quadrant_unmeasured}; it
is not included in \code{keep}. For staged CEM the exemption avoids the entire
L+B campaign; for EM the already fitted banks remain available. The default measures without this heuristic.
The nontrivial-only measurement policy cannot supply trivial anchors and is
therefore incompatible with this option.
}
\section{Output and plotting}{
The result is a \code{CEC_adaptive_grid} with the usual grid components.
\code{summary} contains active leaves only; u0/u1 are log-C edges, v0/v1 lambda
edges, with id, parent and depth. \code{adaptive$mesh$archive} retains parent
geometry. \code{A_harm}, \code{keep_harmonization}, \code{measurement_status}
and \code{quadrant_witness} record additional controls.
For CEM, \code{adaptive$pilot} and \code{adaptive$pilot_fits} preserve the pilot,
\code{measurements} stores campaigns by cell ID, and \code{executions} records
phase, kind, cell, seed, parent-start use and checkpoint reuse. Unmeasured
rows in \code{winners} contain NA. Checkpoints store each pilot fit, each
independent campaign execution and each completed harmonization. Resume
replays the mesh decisions and reuses those numerical results.
Plots draw these actual rectangles, rather than inventing a Cartesian grid.
Extraction at an interior point returns the representative of its containing
cell; points outside the window are rejected. One-C paths describe evaluated
centres at that C, not an interpolated full vertical cut.
}
\seealso{\link{CECfitBoundGrid}, \link{CECrescoreGrid}, \link{CECoptimizers}, \link{CECdiagnostics}.}
\examples{
set.seed(4)
X <- c(rnorm(12,-2,.4), rnorm(12,2,.4))
g <- CECfitBoundGrid(Z=X, grid_strategy="adaptive",
  optimization="em_bank", penalty="bic", k0=2, B=2,
  optimizer_control=list(Rmax=3,Ninit=1,Cmax=10,max_iter=15),
  adaptive_control=list(max_cells=100,max_rounds=2,max_extensions=1))
g$adaptive$stop_reason
CECplotGrid(g)
}
