#' CEClust: Composite Entropy Clustering for Mixed Data
#'
#' `CEClust` provides methods for composite entropy clustering on Gaussian,
#' categorical, and mixed data sets. The package exposes a compact user API for
#' fixed-lambda fits, lambda-grid diagnostics, C by lambda grid exploration, and
#' static or interactive visual inspection.
#'
#' The typical workflow is:
#'
#' 1. configure the runtime with [CECconfigure_runtime()];
#' 2. fit one model with [CECclassif()], scan a lambda grid with
#'    [CECfitLambdaGrid()], or scan a C by lambda grid with
#'    [CECfitBoundGrid()];
#' 3. select stable lambda values with [CECselectStableLambdas()] and extract
#'    representative partitions with [CECextractBestPartitions()];
#' 4. inspect summaries with [CECsummariseGrid()], static plots with
#'    [CECplotGrid()], [CECplotPartition()], and [CECplotPath()], or the Shiny
#'    explorer with [CECexplore()].
#'
#' @section Main entry points:
#' - [CECconfigure_runtime()] configures pure-R or optional fast-backend
#'   execution.
#' - [CECclassif()] fits a repeated-shot CEC model for one `lambda`.
#' - [CECfitLambdaGrid()] evaluates a grid of `lambda` values
#'   using linked forward and backward paths.
#' - [CECfitBoundGrid()] evaluates a C by `lambda` grid.
#' - [CECfitPreset()] runs named preset workflows such as `"iris"`,
#'   `"breast_cancer"`, and `"gaussian_1d"`.
#' - [CECselectStableLambdas()] and [CECsummariseGrid()] provide compact
#'   diagnostics.
#' - [CECplotGrid()], [CECplotPartition()], [CECplotPath()], and [CECexplore()]
#'   provide static and interactive visualisation.
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
