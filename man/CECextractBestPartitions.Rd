% Maintained by hand; see CONTRIBUTING.md.
\name{CECextractBestPartitions}
\alias{CECextractBestPartitions}
\title{Extract representatives along a lambda path}
\description{
Return one fitted representative per lambda from a lambda-diagnostic object. For consolidated grids these representatives have already been selected from the pooled candidates.
}
\usage{
CECextractBestPartitions(lambda_diag, source = c("all", "fixed", "bootstrap"), criterion = c("projected_H", 
    "raw_H"), Z = NULL)
}
\arguments{
\item{lambda_diag}{Output of \code{CECfitLambdaGrid()}, or a grid run component \code{lambda_diag}.}
\item{source}{For solver diagnostics, only \code{"all"} is supported. For historical diagnostics, restrict candidates to original (\code{"fixed"}), bootstrap, or both kinds of runs.}
\item{criterion}{For solver diagnostics, only \code{"projected_H"} is supported: the hard criterion on original observations, including the selected penalty. Historical \code{"raw_H"} ranks training scores and is not equivalent.}
\item{Z}{Original observations for historical projection or decomposition; solver diagnostics already retain the needed fits.}
}
\value{
A list with \code{best} (one representative per lambda), \code{summary}, source/criterion metadata, and for new engines optimizer/penalty metadata.
}
\details{
The result is a representative path, not the independent runs used in
the stability averages. Choosing a representative cannot improve those averages.
Grid compaction can remove repetitions and strip the detailed lambda object.
For a compacted grid, use \code{CECextractPartition()} or the already stored
\code{runs_by_C[[i]]$best_parts}, rather than re-extracting from a stripped
\code{lambda_diag}.
}
\seealso{\link{CECextractPartition}, \link{CECdetectPartitionChanges}, \link{CECresults}.}
\examples{
set.seed(12)
X <- cbind(x = c(rnorm(30, -2, 0.5), rnorm(30, 2, 0.5)),
           y = rnorm(60, sd = 0.5))
d <- CECfitLambdaGrid(X, C = 1, lambda_grid = c(0.6, 1, 1.5),
  optimization = "em_bank", penalty = "none", k0 = 3, B = 3,
  seed = 14, optimizer_control = list(Rmax = 3, Ninit = 1, max_iter = 30),
  verbose = FALSE)
best <- CECextractBestPartitions(d)
best$summary[, c("lambda", "REO", "H_projected")]
}
