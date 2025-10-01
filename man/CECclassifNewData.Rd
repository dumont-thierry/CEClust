\name{CECclassifNewData}
\alias{CECclassifNewData}
\title{Affectation de nouveaux individus à des clusters existants}
\description{
À partir d’un objet \code{params} appris par \code{\link{CECclassif}}, assigne à chaque ligne de \code{Zpred} le cluster le plus probable.
}
\usage{
CECclassifNewData(Zpred, params, idColToPred)
}
\arguments{
  \item{Zpred}{\code{data.frame}/\code{matrix} de nouvelles observations (même schéma que \code{Z} à l’apprentissage, sauf colonnes à prédire).}
  \item{params}{\code{list}. Paramètres du modèle retournés par \code{\link{CECclassif}}.}
  \item{idColToPred}{\code{integer} ou \code{character}. Indices/noms des colonnes ignorées pour l’affectation.}
}
\value{
\code{integer} vector des indices de cluster (longueur \code{nrow(Zpred)}).
}
\seealso{\code{\link{CECclassif}}, \code{\link{CECpredict}}}
\examples{
# suite de l'exemple de \code{CECclassif()}

nnew <- 10
# Deux clusters gaussiens sur 2 variables
x1 <- rnorm(nnew, 0, 1); x2 <- rnorm(nnew, 0, 1)
x3 <- rnorm(nnew, 3, 1); x4 <- rnorm(nnew, 3, 1)
X  <- rbind(cbind(x1, x2), cbind(x3, x4))
# Une variable catégorielle corrélée au cluster
lab <- factor(rep(c("A","B"), each = nnew))

Znew  <- data.frame(X1 = X[,1], X2 = X[,2])
  
pred_clusters <- CECclassifNewData(Znew, fit$params, idColToPred = "Cat")
pred_clusters
}
