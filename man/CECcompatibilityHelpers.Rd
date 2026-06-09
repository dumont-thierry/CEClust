\name{CECcompatibilityHelpers}
\alias{CECaddSmoothedDiagnostics}
\alias{CECdetectPartitionChanges}
\alias{CECextractBestPartitions}
\alias{CECidentifyStableLambdas}
\alias{CECrepairBestPartitionTrajectory}
\title{Compatibility helpers for CEClust grid workflows}
\usage{
CECaddSmoothedDiagnostics(lambda_diag, k = 3,
  ratio_col = "stab_ratio_mean", stab0_col = "stab0_mean",
  stabB_col = "stabB_mean", H_col = "meanH_B")

CECdetectPartitionChanges(best_partitions_obj,
  partition_metric = c("match", "ARI"), sim_threshold = 0.8,
  criterion = c("REO_or_threshold", "REO_only", "threshold_only",
    "REO_and_threshold"),
  include_first = TRUE, method = c("consecutive", "regime_reference"),
  lambda_subset = NULL)

CECextractBestPartitions(lambda_diag, source = c("all", "fixed", "bootstrap"),
  criterion = c("projected_H", "raw_H"), Z = NULL)

CECidentifyStableLambdas(lambda_diag, rule = c("ratio", "stabB", "both"),
  ratio_threshold = 0.8, stabB_threshold = 0.8, sat_threshold = 0.8,
  use_smoothed = TRUE, min_consecutive = 1, best_partitions_obj = NULL,
  ratio_col = "stab_ratio_mean", stabB_col = "stabB_mean",
  sat_col = "sat")

CECrepairBestPartitionTrajectory(best_partitions_obj, Z, tol = 1e-10,
  max_iter = 100)
}
\arguments{
\item{lambda_diag}{Lambda-grid diagnostic object.}
\item{k}{Window size for moving-average smoothing.}
\item{ratio_col, stab0_col, stabB_col, H_col, sat_col}{Column names used in
diagnostic summaries.}
\item{best_partitions_obj}{Object containing best partitions along a lambda
path.}
\item{partition_metric}{Partition-comparison metric.}
\item{sim_threshold}{Similarity threshold used to detect partition changes.}
\item{criterion}{Rule used to detect partition changes.}
\item{include_first}{Whether the first lambda value is marked as a change.}
\item{method}{Change-detection method.}
\item{lambda_subset}{Optional subset of lambda values.}
\item{source}{Which fitted solutions are used to extract best partitions.}
\item{Z}{Reference data.}
\item{rule}{Stability rule.}
\item{ratio_threshold, stabB_threshold, sat_threshold}{Selection thresholds.}
\item{use_smoothed}{Whether smoothed diagnostics are used.}
\item{min_consecutive}{Minimum number of consecutive retained lambdas.}
\item{tol}{Numerical tolerance.}
\item{max_iter}{Maximum number of repair iterations.}
}
\value{
The value returned by the corresponding helper.
}
\description{
These helpers are exported for compatibility with the compact grid engine and
for advanced users who need to inspect lambda-path diagnostics directly. The
recommended entry points remain \code{\link[=CECfitBoundGrid]{CECfitBoundGrid()}},
\code{\link[=CECplotGrid]{CECplotGrid()}}, and
\code{\link[=CECexplore]{CECexplore()}}.
}
