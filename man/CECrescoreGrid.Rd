% Maintained by hand; see CONTRIBUTING.md.
\name{CECrescoreGrid}
\alias{CECrescoreGrid}
\title{Re-evaluate retained EM banks without repeating EM}
\description{
Reuse all independent EM candidate banks on another scoring grid or with another penalty. The winners of each bank, the displayed partitions and the exhaustive diagnostics are recomputed.
}
\usage{
CECrescoreGrid(x,
  C_grid = if (inherits(x, "CEC_adaptive_grid")) x$adaptive$initial_C else x$C_grid,
  lambda_grid = if (inherits(x, "CEC_adaptive_grid")) x$adaptive$initial_lambda else x$lambda_grid,
  penalty = x$penalty, rsi_policy = x$meta$rsi_policy,
  grid_strategy = if (inherits(x, "CEC_adaptive_grid")) "adaptive" else "fixed",
  adaptive_control = if (inherits(x, "CEC_adaptive_grid")) x$adaptive$control else list())
}
\arguments{
\item{x}{A retained \code{CEC_solver_grid} with \code{optimization="em_bank"}. For one-C diagnostics use \code{x$grid}.}
\item{C_grid, lambda_grid}{New distinct positive coordinates; every C must be at most the Cmax used to fit the banks. Defaults to the initial window for adaptive objects, otherwise the original coordinates.}
\item{grid_strategy}{Fixed or adaptive evaluation. Defaults to the strategy of x.}
\item{adaptive_control}{Adaptive controls; defaults to those of an adaptive input. Each penalty can produce its own window and refinement.}
\item{penalty}{\code{"none"} or \code{"bic"}; defaults to the original penalty.}
\item{rsi_policy}{Measurement policy; defaults to the original policy.}
}
\value{
A new \code{CEC_solver_grid} with reused banks and updated summaries, winner indices and diagnostics. The input is unchanged. The result has \code{meta$reused_banks=TRUE}.
}
\details{
No EM descent is repeated. Hard group profiles are still reevaluated
under each C. Bootstrap winners are chosen on their own training data before
comparison on the common original observations. Raw stability keeps these
independent winners. The conditional reference rule is checked before filtering
each measured point and can run new single-descent CEM bootstrap fits, without
repeating EM. Its setting is inherited from the input. See \link{CECdiagnostics}.

Stored data are checked against their original fingerprint. Do not modify Z,
bank contents or bootstrap indices before rescoring. The original randomization,
measurement-band limit and filtering thresholds are preserved. Use the selection
and plotting functions for other display thresholds, or a new fitting call with
a matching checkpoint to change the measurement-band limit.
This operation cannot add candidates absent from the fitted banks or extend C
beyond Cmax. It does not save a result file automatically.
}
\seealso{\link{CECfitBoundGrid}, \link{CECoptimizers}, \link{CECresults}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
g <- CECfitBoundGrid(Z = X, C_grid = c(0.2, 1),
  lambda_grid = c(0.6, 1, 1.5), optimization = "em_bank",
  penalty = "none", k0 = 3, B = 3, algo_seed = 14,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  n_cores = 1)
g_bic <- CECrescoreGrid(g, penalty = "bic")
cbind(CEC = g$summary$REO, CEC_BIC = g_bic$summary$REO)
stopifnot(identical(g$banks, g_bic$banks))
}
