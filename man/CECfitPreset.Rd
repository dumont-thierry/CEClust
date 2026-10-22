% Maintained by hand; see CONTRIBUTING.md.
\name{CECfitPreset}
\alias{CECfitPreset}
\title{Fit a named demonstration dataset}
\description{
Convenience entry points for a generated one-dimensional mixture, iris, and the BreastCancer data in mlbench. For a controlled analysis of other data, use CECfitBoundGrid directly.
}
\usage{
CECfitPreset(preset = c("gaussian_1d", "iris", "breast_cancer"), ..., save_results = FALSE, 
    output_dir = file.path(tempdir(), "CEClust"), checkpoint_dir = FALSE, auto_checkpoint = FALSE, 
    resume = FALSE, force_recompute = FALSE)
}
\arguments{
\item{preset}{\code{"gaussian_1d"}, \code{"iris"}, or \code{"breast_cancer"}.}
\item{...}{Grid arguments, including explicit optimizer, penalty and budgets. Iris and breast-cancer representation options are described below.}
\item{save_results}{Logical; save the result as an RDS file. The default is \code{FALSE}.}
\item{output_dir}{Directory for requested results or automatic checkpoints. The default temporary directory is not persistent between R sessions.}
\item{checkpoint_dir}{\code{FALSE} disables checkpoints. A directory path enables them; \code{TRUE} uses a subdirectory of \code{output_dir}. See the checkpoint section below.}
\item{auto_checkpoint}{If \code{TRUE} and no directory is supplied, use \code{file.path(output_dir,"solver_checkpoints")} for the consolidated engine.}
\item{resume}{If \code{TRUE}, reuse a complete matching fitting checkpoint when present. The default here is \code{FALSE}.}
\item{force_recompute}{If \code{TRUE}, ignore an existing fitting checkpoint.}
}
\value{
A \code{CEC_bound_grid}; explicit new-engine calls also inherit from \code{CEC_solver_grid}. Display labels and preprocessing metadata are retained.
}
\details{
The iris preset excludes species from fitting unless
\code{include_species_as_feature=TRUE}. The breast-cancer preset excludes the
known class unless \code{include_class_as_feature=TRUE}. Adding a factor feature
requires the historical mixed-data CEM, not a new numeric engine.

\code{quantitative_representation} accepts \code{"raw"}, \code{"normalized"},
\code{"PC2"} through the available dimension, and corresponding \code{"rawPC2"}
forms. PC representations use standardized variables; rawPC forms do not.
Iris defaults to raw measurements, breast cancer to normalized measurements.
The latter requires \pkg{mlbench} and prepares complete quantitative rows.

Preset-specific defaults can override full-grid defaults: iris and breast cancer
use 20 original and 20 bootstrap repetitions, 50 fresh starts, eight workers
and thresholds 0.8, 0.8, 0.05. State budgets and thresholds explicitly when
comparing datasets or optimizers.
}
\seealso{\link{CECfitBoundGrid}, \link{CECoptimizers}.}
\examples{
g <- CECfitPreset("iris", C_grid = 1, lambda_grid = c(0.8, 1.2),
  optimization = "em_bank", penalty = "none", k0 = 2, B = 1,
  optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 20),
  algo_seed = 14, n_cores = 1, verbose = FALSE,
  stab_algo_threshold = 0.85, rsi_threshold = 0.95, sat_threshold = 0)
CECsummariseGrid(g)
}
