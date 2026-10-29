\name{CECdiagram}
\alias{CECdiagram}
\alias{CECdiagramControl}
\alias{CECinitialGrid}
\title{Choose Between Adaptive Grids and Candidate Bands}
\usage{
CECdiagram(X, method = c("adaptive", "bands"), L = 20L, B = 20L,
  workers = 1L, seed = 1L, penalty = c("none", "bic"),
  control = list(), C_grid = NULL, lambda_grid = NULL,
  bootstrap_indices = NULL, directory = tempfile("CEC-diagram-"),
  resume = FALSE)
CECdiagramControl(method = c("adaptive", "bands"), ...)
CECinitialGrid(X, powers = 1:10)
}
\arguments{
\item{X}{Finite numeric matrix or numeric data frame, observations in rows.}
\item{method}{Search engine. Adaptive grid is the default; bands remains experimental.}
\item{L, B}{Numbers of independent original-data and bootstrap repetitions; at least two and one, respectively.}
\item{workers}{Number of worker processes.}
\item{seed}{Nonnegative integer seed.}
\item{penalty}{Unpenalized or Gaussian BIC-type order penalty.}
\item{control}{Named list of settings from CECdiagramControl.}
\item{C_grid, lambda_grid}{Optional initial centres for adaptive grids. C centres must be log-equally spaced and lambda centres linearly spaced. Not accepted for bands.}
\item{bootstrap_indices}{Optional list of B bootstrap index vectors.}
\item{directory}{Archive directory for checkpoints and diagram.rds. Preserve it for the bands result, which references files on disk.}
\item{resume}{Resume compatible checkpoints. Historical standalone-package archives have different runtime signatures and are not reused.}
\item{...}{Method-specific settings. Adaptive accepts arguments of the underlying adaptive engine, such as optimizer_control, adaptive_control, optimization, rsi_threshold and reference_bootstrap. Bands accepts arguments of bands_control. Unknown names are errors.}
\item{powers}{Increasing equally spaced integers used as powers of two.}
}
\description{Runs one of two separately implemented search engines within one package. Original lower-level interfaces remain available.}
\details{
For adaptive grids, the default optimization is CEM. Unless C_grid is supplied,
CECinitialGrid computes the Gaussian density peak Cmin from the maximum-likelihood
covariance (divisor n) of the single-group fit. Initial centres are
\eqn{2^k C_{\min}} for k from 1 to 10.
A singular covariance or unrepresentable peak produces an explicit error,
not a hidden regularization.
Cmin is an exploration reference, not a proof of saturation for every partition.
Left extension is disabled by the common interface. The existing rectangular-cell
geometry is retained, including half-cell edges around the initial centres.
Right and both lambda extensions retain their existing rules and finite budgets.
Observed left-boundary saturation remains recorded but is not a required stopping
condition. To restore historical left extension, specify
\code{adaptive_control = list(extend_left = TRUE)} in control.

The historical CECfitBoundGrid interface retains its earlier default initial grid
and four-sided exploration for compatibility. Use CECdiagram for the new initialization.

Bands uses the latest agreed search budget: orders 1:20, five starts, two extra
starts, ten anchors, enriched order at most 40, 11 initial and 41 maximum C columns,
two refinement and two projection rounds, 200 projection queries. Other bands
defaults, including the strict RSI threshold 0.95, are preserved.
Direct cec_bands and bands_control calls retain their historical defaults.
No automatic cleanup is added: bands_cleanup is a separate removal-only audit.

BIC is an order penalty within CEC, not a separate density-model selector.
The engines retain their different stability and reference-selection procedures.
Switching engines is not guaranteed to give the same partitions.
}
\value{
CECdiagram returns the native result of the chosen engine, with an additional
CECdiagram attribute recording its method, version, archive and automatic
initialization. Adaptive results use CECsummariseGrid, CECextractBestPartitions
and CECplotGrid. Bands results use bands_partitions, bands_plot and bands_validate.
CECdiagramControl returns a named list. CECinitialGrid returns Cmin, C_grid and powers.
}
\examples{
X <- matrix(seq(-2, 2, length.out = 40), ncol = 1)
CECinitialGrid(X)
CECdiagramControl("bands", max_anchors = 5L)
}
\seealso{\code{\link{CECfitBoundGrid}}, \code{\link{cec_bands}}, \code{\link{bands_cleanup}}}
