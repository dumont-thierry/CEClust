% Maintained by hand; see CONTRIBUTING.md.
\name{CECselectStableLambdas}
\alias{CECselectStableLambdas}
\title{Filter one C slice by stability and saturation}
\description{
Identify all admissible lambda cells at a selected C. The function returns a filtered view; it does not impose a unique clustering or change the fitted grid.
}
\usage{
CECselectStableLambdas(x, best_parts = NULL, C = NULL, ...)
}
\arguments{
\item{x}{A full grid, a one-C grid result, or lambda diagnostics. Consolidated lambda diagnostics must retain their parent \code{grid}.}
\item{best_parts}{Historical best-partitions object, used to recover saturation when x contains only legacy diagnostics. Not needed for a solver grid.}
\item{C}{C slice; nearest stored C if not exact. If omitted, the display-default cell determines C; see details.}
\item{...}{\code{stab_algo_threshold}, \code{rsi_threshold}, \code{sat_threshold}, \code{use_smoothed}, \code{min_consecutive}. For consolidated grids, NULL thresholds inherit metadata, smoothing must be FALSE, and min_consecutive defaults to 1.}
}
\value{
A list with \code{summary}, \code{stable_lambda}, \code{stable_index}, logical \code{stable}, rejection information and the applied thresholds. \code{lambda_min} is the first retained lambda or NA; it is not an automatically selected analysis scale.
}
\details{
For the consolidated engine, retention requires available finite
diagnostics, Stab_algo strictly above its threshold, RSI strictly above its
threshold, and saturation at or below its threshold. In particular,
sat_threshold=0 allows exactly zero saturation. Missing diagnostics never pass.
Changing a threshold cannot recover a diagnostic skipped at fitting time;
for EM banks, use \code{CECrescoreGrid()} with a broader measurement policy.
Adaptive component diagrams can use values inherited from a passing witness;
this is not a measurement at every retained cell. Inspect \code{diagnostic_origin}
and \code{diagnostic_witness}. Changed thresholds may require a new traversal;
see \link{CECadaptiveDiagram}.

\code{min_consecutive} removes retained runs shorter than this many adjacent
grid entries. It counts cells, not a distance in lambda. If coordinates are omitted, the display selects a retained cell with smallest REO, taking the first tie. If no cell is retained, it uses a finite-REO cell instead. This is a display default, not a rule for choosing a final clustering.

Historical selectors have their own inclusive stability comparisons and defaults.
They are documented in \link{CECmigration}; do not use their output as if it
were the same rule as the consolidated selector.
}
\seealso{\link{CECdiagnostics}, \link{CECrescoreGrid}, \link{CECplotPath}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
g <- CECfitBoundGrid(Z = X, C_grid = c(0.2, 1),
  lambda_grid = c(0.6, 1, 1.5), optimization = "em_bank",
  penalty = "none", k0 = 3, B = 3, algo_seed = 14,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  n_cores = 1)
selected <- CECselectStableLambdas(g, C = 1,
  stab_algo_threshold = 0.85, rsi_threshold = 0.95, sat_threshold = 0)
selected$stable_lambda
}
