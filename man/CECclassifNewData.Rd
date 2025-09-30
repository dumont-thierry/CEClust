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
Znew <- Z[1:5, ]
Znew$Cat <- factor(NA, levels = levels(Z$Cat))
pred_clusters <- CECclassifNewData(Znew, fit$params, idColToPred = "Cat")
pred_clusters
}
