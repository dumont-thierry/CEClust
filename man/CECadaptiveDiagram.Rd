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
\item{\code{harmonize=TRUE, harmonization_min=0.95}}{For CEM, propagate improving neighbour partitions and run a descent from them. The initial pass visits the whole mesh. Later passes start from the newly evaluated cells and their neighbours; an improvement reactivates affected neighbours until propagation stops. Keep requires agreement at least harmonization_min between the pre-harmonization cell winner and its final representative. Independent repetitions remain unchanged. EM uses its common bank, has no additional harmonization, and has A_harm=1.}
\item{\code{trivial_quadrants=FALSE}}{Optional heuristic exemption from diagnostic measurement in validated upper-left trivial quadrants.}
\item{\code{audit_fraction=0.10, anchor_stability=0.95}}{Quadrant controls described below.}
}
Inspect \code{g$adaptive$history}, \code{margins}, \code{restrictions},
\code{stop_reason} and \code{completed}. A budget stop or a Cmax stop is never
reported as completed adaptation. Even completed adaptation only means that
the specified numerical stopping rules were met within the explored window.
}
\section{Diagnostics and optional quadrants}{
Measurement policies and exhaustive comparisons are unchanged. EM fits are
shared across cells within each independent execution, but winners remain
independent between executions. With CEM, each evaluated cell has its own
independent repetitions. Neighbour harmonization never substitutes its output
into those repetitions. Final diagnostics use the final representative's
saturation to determine the measurement policy.
Identical winner patterns and partition comparisons are cached exactly.

Optional trivial quadrants use a measured, trivial anchor passing RSI and
algorithmic-stability thresholds. Every explored cell with C no greater and
lambda no smaller must be trivial. Controls use one deterministic, persistent
budget for the union of candidate quadrants; overlapping quadrants do not
duplicate that budget. A failing control invalidates its quadrant. An exempted
cell retains missing diagnostics and is marked \code{quadrant_unmeasured}; it
is not included in \code{keep}. The default measures without this heuristic.
The nontrivial-only measurement policy cannot supply trivial anchors and is
therefore incompatible with this option.
}
\section{Output and plotting}{
The result is a \code{CEC_adaptive_grid} with the usual grid components.
\code{summary} contains active leaves only; u0/u1 are log-C edges, v0/v1 lambda
edges, with id, parent and depth. \code{adaptive$mesh$archive} retains parent
geometry. \code{A_harm}, \code{keep_harmonization}, \code{measurement_status}
and \code{quadrant_witness} record additional controls.
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
