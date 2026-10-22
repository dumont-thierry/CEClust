% Maintained by hand; see CONTRIBUTING.md.
\name{CECfitLambdaGrid}
\alias{CECfitLambdaGrid}
\title{Fit and diagnose a lambda path at one C}
\description{
One-C interface to independent CEM searches or EM candidate banks. It returns lambda diagnostics and their representative partitions. Historical linked trajectories remain available for compatibility.
}
\usage{
CECfitLambdaGrid(Z, lambda_grid = seq(0.1, 2, by = 0.1), k0 = 10, B = 10, C = 1, Cquali = Inf, 
    r0 = NULL, Nshots_fresh = 5, Nshots_warm = 5, Nshots_init = NULL, Nloop = 300, familyType = "gaussAndDiscreteVector", 
    sizeMaxOutlier = 0, autoRegroupOutliers = FALSE, focus = NULL, replace = TRUE, seed = NULL, 
    silent = TRUE, verbose = TRUE, n_cores = 1L, checkpoint_dir = NULL, auto_checkpoint = NULL, 
    resume = TRUE, batch_size = NULL, show_progress = NULL, rerun_failed_serial = TRUE, 
    checkpoint_compress = FALSE, partition_metric = c("match", "ARI"), stability_approach = c("exhaustive", 
        "pairwise", "reference_bestH"), optimization = c("cem", "em_bank"), penalty = c("none", 
        "bic"), optimizer_control = list(), rsi_policy = c("all", "unsaturated", "nontrivial_unsaturated", 
        "saturation_band"), rsi_saturation_max = 0.1)
}
\arguments{
\item{Z}{Data as for \code{CECclassif()}.}
\item{lambda_grid}{Distinct positive lambda values, sorted for the consolidated engine.}
\item{k0}{Number L of original repetitions; default 10.}
\item{B}{Number of bootstrap repetitions; default 10.}
\item{C}{Positive density ceiling.}
\item{Cquali}{Historical qualitative bound; new engines require Inf.}
\item{r0}{Default Rmax for numeric engines, or historical initial group count. If NULL, the new engine uses min(20,n).}
\item{Nshots_fresh}{Default Ninit for numeric engines; fresh starts per historical trajectory step.}
\item{Nshots_warm}{Historical warm starts; ignored by the consolidated engine.}
\item{Nshots_init}{Compatibility alias used when Nshots_fresh is omitted.}
\item{Nloop}{Default max_iter for numeric engines; historical descent limit.}
\item{familyType}{Reference family; see \link{CECclassif}.}
\item{sizeMaxOutlier, autoRegroupOutliers, focus}{Historical controls; use 0, FALSE and NULL with new engines.}
\item{replace}{Bootstrap sampling with replacement. Must be TRUE for new engines.}
\item{seed}{Master seed for this interface. If NULL, the consolidated route draws a seed from the caller RNG. A control seed takes precedence.}
\item{silent}{Historical fit-output switch; not used by the consolidated engine.}
\item{verbose}{Print workflow information.}
\item{n_cores}{Number of local workers for independent repetitions.}
\item{checkpoint_dir}{FALSE/NULL disables consolidated checkpoints, a path enables them, TRUE chooses a temporary directory. Historical behavior can differ; set it explicitly for reproducibility.}
\item{auto_checkpoint}{Explicit TRUE enables consolidated temporary checkpoints if no path is supplied. In the historical route NULL enables them for parallel runs.}
\item{resume}{Reuse matching completed starts and repetitions, or a complete fitting cache; default TRUE for this interface.}
\item{batch_size, show_progress, rerun_failed_serial, checkpoint_compress}{Historical task scheduling, progress, retry and compression controls. Not used by the consolidated engine.}
\item{partition_metric}{Optimal label matching or adjusted Rand index; see \link{CECstability}.}
\item{stability_approach}{Default \code{"exhaustive"}: original-original and original-bootstrap pairs. Historical \code{"pairwise"} uses bootstrap-bootstrap pairs for bootstrap stability; \code{"reference_bestH"} uses a selected reference. New engines require exhaustive.}
\item{optimization, penalty, optimizer_control}{As in \link{CECclassif}. Explicit choices select the consolidated engine.}
\item{rsi_policy, rsi_saturation_max}{Measurement policy and band limit as in \link{CECfitBoundGrid}.}
}
\value{
The consolidated route returns \code{CEC_solver_lambda}: \code{summary}
uses columns \code{lambda}, \code{stab0_mean}, \code{stabB_mean},
\code{stab_ratio_mean}, \code{sat} and \code{measured}. It also contains
\code{best_parts} and the full parent \code{grid}. Use \code{d$grid} with
grid plotting or rescoring functions, and \code{CECextractBestPartitions(d)}
for a path of representatives. Historical results store per-lambda repetitions
in \code{details}.
}
\details{
The numeric route uses the same computations as \code{CECfitBoundGrid()}
at one C. Its initial selection thresholds are 0.85, 0.95 and zero; change the
view afterwards through \code{CECselectStableLambdas()} or use the full-grid
interface to specify thresholds at fitting time. Comparisons remain exhaustive
and unsmoothed. Specifying only \code{rsi_saturation_max} does not switch engines;
state \code{rsi_policy} explicitly.
}
\seealso{\link{CECfitBoundGrid}, \link{CECselectStableLambdas}, \link{CECmigration}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
d <- CECfitLambdaGrid(X, C = 1, lambda_grid = c(0.6, 1, 1.5),
  optimization = "em_bank", penalty = "none", k0 = 3, B = 3,
  seed = 14, optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  verbose = FALSE)
d$summary
best <- CECextractBestPartitions(d)
best$summary[, c("lambda", "REO", "H_projected")]
CECplotGrid(d$grid)
}
