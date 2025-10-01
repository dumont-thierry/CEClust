#' CEClust: Composite Entropy Clustering for mixed data
#'
#' @description
#' CEClust provides tools for unsupervised clustering of mixed-type data
#' (numeric and categorical) using the Composite Entropy Criterion (CEC).
#'
#' The method fits mixtures with Gaussian components for numeric variables and
#' multinomial components for categorical variables, optimized under a composite
#' entropy objective. Main user functions:
#'
#' \itemize{
#'   \item \code{\link{CECclassif}} — learn a clustering model from data;
#'   \item \code{\link{CECclassifNewData}} — assign new observations to existing clusters;
#'   \item \code{\link{CECpredict}} — predict missing columns conditionally on the cluster.
#' }
#'
#' @section Dependencies:
#' Imports \pkg{stats} and \pkg{mvtnorm}.
#'
#' @name CEClust
#' @aliases CEClust-package
#' @docType package
#'
#' @importFrom stats dnorm runif rnorm cov2cor
#' @importFrom mvtnorm dmvnorm
"_PACKAGE"
