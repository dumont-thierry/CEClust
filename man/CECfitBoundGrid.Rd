% Maintained by hand; see CONTRIBUTING.md.
\name{CECfitBoundGrid}
\alias{CECfitBoundGrid}
\title{Fit and diagnose a C by lambda grid}
\description{
Evaluate clustering across a supplied grid, repeat the algorithm independently on the original data and bootstrap samples, and retain cells that meet the requested stability and saturation thresholds.
}
\usage{
CECfitBoundGrid(..., save_results = FALSE, output_dir = file.path(tempdir(), "CEClust"), checkpoint_dir = FALSE, 
    auto_checkpoint = FALSE, resume = FALSE, force_recompute = FALSE, optimization = c("cem", 
        "em_bank"), penalty = c("none", "bic"), optimizer_control = list(), rsi_policy = c("all", 
        "unsaturated", "nontrivial_unsaturated", "saturation_band"), rsi_saturation_max = 0.1, grid_strategy = c("fixed", "adaptive"),
    adaptive_control = list(), reference_bootstrap = FALSE,
    require_reference_agreement = NULL)
}
\arguments{
\item{reference_bootstrap}{Opt into the conditional fixed-reference bootstrap rule. The default FALSE preserves independent-only diagnostics. See \link{CECdiagnostics}.}
\item{require_reference_agreement}{NULL enables the independent-reference gate unless conditional bootstraps are explicitly requested. TRUE requires agreement at least stab_algo_threshold between the best independent original partition and the grid reference. FALSE disables this gate. Adaptive controls override these settings. Explicit TRUE cannot be combined with reference_bootstrap=TRUE.}
\item{...}{Named data, grid, fitting and diagnostic arguments listed below. For the recommended workflow, supply \code{optimization} explicitly. Unimplemented options are rejected by the consolidated engine.}
\item{save_results}{Logical; save the result as an RDS file. The default is \code{FALSE}.}
\item{output_dir}{Directory for requested results or automatic checkpoints. The default temporary directory is not persistent between R sessions.}
\item{checkpoint_dir}{\code{FALSE} disables checkpoints. A directory path enables them; \code{TRUE} uses a subdirectory of \code{output_dir}. See the checkpoint section below.}
\item{auto_checkpoint}{If \code{TRUE} and no directory is supplied, use \code{file.path(output_dir,"solver_checkpoints")} for the consolidated engine.}
\item{resume}{If \code{TRUE}, reuse matching completed starts and repetitions, or a complete fitting cache when present. The default here is \code{FALSE}.}
\item{force_recompute}{If \code{TRUE}, ignore an existing fitting checkpoint.}
\item{optimization}{\code{"cem"} or \code{"em_bank"}. Explicitly supplying this argument selects the consolidated numeric Gaussian engine.}
\item{penalty}{\code{"none"} or \code{"bic"}. Explicitly supplying it also selects the consolidated engine.}
\item{optimizer_control}{Named controls described in \link{CECoptimizers}. A nonempty list selects the consolidated engine.}
\item{rsi_policy}{\code{"all"}, \code{"unsaturated"}, \code{"nontrivial_unsaturated"}, or \code{"saturation_band"}. Controls where both stability measures and RSI are computed. Adaptive CEM defaults to \code{"saturation_band"} when omitted and applies the policy before launching cell campaigns; other routes default to \code{"all"}. See \link{CECdiagnostics}. Explicitly supplying it selects the consolidated engine.}
\item{grid_strategy}{\code{"fixed"} (default) or \code{"adaptive"}. Adaptive mode integrates four-sided extension and frontier refinement; see \link{CECadaptiveDiagram}.}
\item{adaptive_control}{Named adaptive-window, refinement, harmonization and diagnostic controls. Adaptive diagrams default to \code{diagnostic_strategy="components"}; use \code{"cells"} for cell-by-cell coverage. See \link{CECadaptiveDiagram}.}
\item{rsi_saturation_max}{Maximum saturation for the measurement policy \code{"saturation_band"}; default 0.10. It is separate from the final \code{sat_threshold}. Used by adaptive CEM's default policy; elsewhere set \code{rsi_policy="saturation_band"} explicitly.}
}
\value{
A \code{CEC_solver_grid} inheriting from \code{CEC_bound_grid} and
\code{CEC_grid}, with \code{summary}, \code{meta}, data, candidate banks,
independent winners and representative partitions. Calls using the historical
engine return \code{CEC_bound_grid} without the solver-specific components.
See \link{CECresults} for the object structure.
}
\section{Data and grid arguments}{
The following named arguments are accepted through \code{...} by the consolidated engine.
\describe{
\item{\code{Z}}{Finite complete numeric data, observations in rows. No automatic standardization. If omitted, a one-dimensional Gaussian-mixture demonstration is generated, using \code{n=100}, \code{data_seed=13}.}
\item{\code{C_grid, lambda_grid}}{Distinct, finite, strictly positive coordinates, sorted internally. Defaults are \code{1} and \code{seq(0.5,2,by=0.1)}. In fixed mode these are evaluated coordinates. In adaptive mode they are the initial centres; omit them for a data-scaled initial window, or supply log-regular C and regular lambda coordinates. See \link{CECadaptiveDiagram}.}
\item{\code{display_data, labels, label_name, dataset_name, density}}{Optional display information. Rows must match Z. External labels are not used to fit the clustering. A supplied density is for a one-dimensional display.}
\item{\code{familyType, Cquali}}{Optional family metadata; new engines require numeric Gaussian data and \code{Cquali=Inf}.}
}
}
\section{Budgets, randomization and diagnostics}{
\describe{
\item{\code{k0, B}}{Original and bootstrap repetitions, corresponding to L and B in the diagnostic definitions. Both default to 20; k0 must be positive, B can be zero. At least two original repetitions and one bootstrap repetition are required for a finite RSI.}
\item{\code{r0, Nshots_fresh, Nloop}}{Fallbacks for \code{Rmax}, \code{Ninit} and \code{max_iter}; defaults 20, 5 and 1000. Explicit \code{optimizer_control} fields take precedence.}
\item{\code{algo_seed}}{Nonnegative integer, default 1. Original runs, bootstrap indices and bootstrap fits receive derived seeds. A non-NULL \code{optimizer_control$seed} takes precedence.}
\item{\code{n_cores}}{Number of local workers, default 1. Independent repetitions are parallelized. With fixed seeds, one and multiple workers use the same randomizations.}
\item{\code{partition_metric}}{\code{"match"} (default) or \code{"ARI"}; see \link{CECstability}.}
\item{\code{stab_algo_threshold, rsi_threshold, sat_threshold}}{Defaults 0.85, 0.95 and 0. Retention requires strict \code{>} for both stability thresholds and \code{<=} for saturation.}
\item{\code{bootstrap_indices}}{Optional list of B integer vectors of length n, with entries in 1:n. Useful for paired comparisons using exactly the same resamples.}
\item{\code{verbose}}{Print the number of cells, repetitions and retained cells; default FALSE.}
}
}
\section{Independent repetitions and the displayed partition}{
Each repetition selects its own winner. Bootstrap winners are selected on
their training sample and projected onto the original observations with frozen
learned parameters. Exhaustive diagnostics compare these independent outputs.
The final displayed partition is selected separately: for EM banks, from all
projected candidates pooled across banks; for CEM, from the independent winners
at that cell on fixed grids. Adaptive CEM defaults to a frozen harmonized
pilot and component traversal. With \code{diagnostic_strategy="cells"}, better
campaign results can update its representative. Pooled or harmonized winners
are never substituted into the raw independent diagnostic repetitions.
The conditional reference branch is recorded separately; see \link{CECdiagnostics}.

EM fits are shared across all cells within one execution; hard profiles are
reevaluated under each C. Fixed-grid CEM searches run separately at each cell.
Adaptive CEM constructs and harmonizes a pilot before running independent
campaigns at visited eligible points. The default component traversal stops
at its first passing witness and inherits its diagnostic values within that
region for display and selection; direct measurements remain identifiable.
On fixed grids neither route applies a final cross-cell repair. EM uses its
common candidate bank and needs no CEM harmonization.
}
\section{Checkpoints and saved results}{
Checkpoints record each completed independent repetition and, for EM, each completed start with its random-generator state. Interrupted fitting can reuse these partial records. A complete bank cache is also written after all repetitions finish. A matching complete checkpoint avoids repeating the fitting stage. Adaptive CEM first constructs its pilot without L+B campaigns. It checkpoints
each pilot fit, each original or bootstrap campaign execution, and each
completed harmonization separately. Mesh decisions are replayed on resume;
completed numerical fits and propagation are reused. Independent campaign
seeds depend on cell ID, kind and repetition, not worker count.
Its identity includes the package version, data, bootstrap indices, budgets
and seeds. For EM banks, it does not depend on the scoring grid or penalty;
for CEM it does. Diagnostics and final scores are computed again on reuse.

Use a persistent \code{checkpoint_dir} for reuse across sessions. A saved result
contains the complete returned object and can be read with \code{readRDS()}.
Automatic result names identify the fitting cache; rescoring the same EM cache
on a new grid can reuse that name. Use explicit distinct \code{saveRDS()} names
to preserve several scoring configurations.
}
\section{Compatibility}{
When no optimization, penalty, RSI-policy or grid-strategy argument is
explicitly supplied and both control lists are empty, the historical
linked-path engine is used.
Its defaults and repair behavior differ; see \link{CECmigration}.
For a new analysis, always state the optimizer and penalty.

\code{Nshots_warm}, batching/progress switches and legacy compactness switches
are accepted but do not alter fitting in the consolidated engine. Warm starts
are not shared between independent repetitions. Setting
\code{final_descent_from_best=TRUE} or \code{use_smoothed_keep=TRUE} is rejected.
Use \code{grid_strategy="adaptive"} for the integrated extension and refinement workflow, including CEM or EM banks and either penalty. See \link{CECadaptiveDiagram}.
}
\seealso{\link{CECclassif}, \link{CECdiagnostics}, \link{CECresults}, \link{CECrescoreGrid}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
g <- CECfitBoundGrid(Z = X, C_grid = c(0.2, 1),
  lambda_grid = c(0.6, 1, 1.5), optimization = "em_bank",
  penalty = "none", k0 = 3, B = 3, algo_seed = 14,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  n_cores = 1)
s <- CECsummariseGrid(g)
s[, c("C", "lambda", "REO", "sat", "stab_algo", "stab_boot", "rsi", "keep")]
CECplotGrid(g)
}
