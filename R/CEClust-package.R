#' CEClust: Composite Entropy Clustering for mixed data
#'
#' @description
#' The CEClust package provides tools for unsupervised classification
#' of mixed-type data (numeric and categorical) based on the
#' Composite Entropy Criterion (CEC).
#'
#' The main exported functions are:
#' \itemize{
#'   \item \code{\link{CECclassif}}: learn a clustering model from data.
#'   \item \code{\link{CECclassifNewData}}: assign new observations
#'         to previously estimated clusters.
#'   \item \code{\link{CECpredict}}: predict missing columns conditionally
#'         on the cluster assignment.
#' }
#'
#' @details
#' CEClust implements Composite Entropy Clustering, where clusters are
#' estimated by minimizing an entropic criterion. The method can handle
#' both continuous (Gaussian) and categorical (multinomial) variables.
#'
#' @docType package
#' @name CEClust
#' @keywords package
#'
#' @importFrom stats dnorm runif rnorm cov2cor
#' @importFrom mvtnorm dmvnorm
NULL
