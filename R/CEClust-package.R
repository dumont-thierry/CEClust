#' CEClust: Composite Entropy Clustering for Mixed Data
#'
#' `CEClust` provides methods for composite entropy clustering on Gaussian,
#' categorical, and mixed data sets. The package covers both repeated-shot
#' estimation for a fixed regularization value and linked-lambda diagnostics to
#' study how cluster structure evolves over a grid of `lambda` values.
#'
#' The typical workflow is:
#'
#' 1. configure the runtime with [CECconfigure_runtime()];
#' 2. fit one model with [CECclassif()] or diagnose a grid with
#'    [CECdiagnose_lambda_grid_linked()];
#' 3. extract representative partitions with [CECextractBestPartitions()];
#' 4. inspect stability and visual summaries with
#'    [CECidentifyStableLambdas()], [plotCECdiagnoseLambdaGrid()], and
#'    [plotCECdiagnoseInteractive()].
#'
#' @section Main entry points:
#' - [CECconfigure_runtime()] configures pure-R or optional fast-backend
#'   execution.
#' - [CECclassif()] fits a repeated-shot CEC model for one `lambda`.
#' - [CECdiagnose_lambda_grid_linked()] evaluates a grid of `lambda` values
#'   using linked forward and backward paths.
#' - [CECextractBestPartitions()] summarises one representative partition per
#'   `lambda`.
#' - [plotCECdiagnoseInteractive()] launches an interactive partition explorer.
#'
#' @aliases CEClust
#' @keywords internal
#' @import graphics
#' @import grDevices
#' @import stats
#' @import utils
#' @useDynLib CEClust, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom mvtnorm dmvnorm
"_PACKAGE"
