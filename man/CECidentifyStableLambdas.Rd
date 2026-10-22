% Maintained by hand; see CONTRIBUTING.md.
\name{CECidentifyStableLambdas}
\alias{CECidentifyStableLambdas}
\title{Apply the historical lambda-selection rule}
\description{
Compatibility selector for historical summaries. New analyses should use CECselectStableLambdas, which also checks algorithmic stability and honors the consolidated grid metadata.
}
\usage{
CECidentifyStableLambdas(lambda_diag, rule = c("ratio", "stabB", "both"), ratio_threshold = 0.8, stabB_threshold = 0.8, 
    sat_threshold = 0.8, use_smoothed = TRUE, min_consecutive = 1, best_partitions_obj = NULL, 
    ratio_col = "stab_ratio_mean", stabB_col = "stabB_mean", sat_col = "sat")
}
\arguments{
\item{lambda_diag}{Historical lambda diagnostics.}
\item{rule}{\code{"ratio"}, \code{"stabB"} or \code{"both"}; compares RSI, bootstrap stability, or both.}
\item{ratio_threshold, stabB_threshold}{Inclusive lower thresholds for the chosen historical diagnostic rule.}
\item{sat_threshold}{Strict upper saturation bound in this historical helper. A value of zero rejects zero saturation too; use the maintained selector for a zero-saturation rule.}
\item{use_smoothed}{Use columns ending in _smooth, which must already exist. Default TRUE is historical.}
\item{min_consecutive}{Minimum length of retained runs in adjacent entries.}
\item{best_partitions_obj}{Optional representatives from which to calculate a missing saturation column.}
\item{ratio_col, stabB_col, sat_col}{Names of the corresponding summary columns.}
}
\value{
A list of retained lambdas, masks, runs and a diagnostic summary.
}
\details{
This historical helper does not impose a separate Stab_algo threshold.
Missing saturation is allowed by its old rule. Its definitions and default
smoothing therefore differ from the consolidated selector; it should not be
used to reproduce the current diagram filter.
}
\seealso{\link{CECselectStableLambdas}, \link{CECmigration}.}
\examples{
old <- list(summary = data.frame(lambda = c(0.6, 1, 1.4),
  stab_ratio_mean = c(0.7, 0.96, 1), stabB_mean = c(0.6, 0.9, 1),
  sat = c(0.2, 0, 0), n_sat_clusters = c(1, 0, 0),
  any_bound_saturation = c(TRUE, FALSE, FALSE)))
CECidentifyStableLambdas(old, use_smoothed = FALSE,
  ratio_threshold = 0.95, sat_threshold = 0.05)$stable_lambda
}
