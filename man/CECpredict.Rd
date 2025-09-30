\name{CECpredict}
\alias{CECpredict}
\title{Prédiction de colonnes manquantes conditionnellement au cluster}
\description{
En utilisant les paramètres \code{params} appris par \code{\link{CECclassif}}, prédit les valeurs des colonnes \code{idColToPred} pour \code{Zpred}.
Les colonnes factor sont prédites par l’argmax des probabilités discrètes; les colonnes numériques via la moyenne conditionnelle gaussienne.
}
\usage{
CECpredict(Zpred, params, idColToPred)
}
\arguments{
  \item{Zpred}{\code{data.frame}/\code{matrix} de nouvelles observations.}
  \item{params}{\code{list}. Paramètres du modèle retournés par \code{\link{CECclassif}}.}
  \item{idColToPred}{\code{integer} ou \code{character}. Colonnes à prédire.}
}
\value{
\code{data.frame} de taille \code{nrow(Zpred) \u00D7 length(idColToPred)} avec les valeurs prédites (facteurs et/ou numériques).
}
\seealso{\code{\link{CECclassif}}, \code{\link{CECclassifNewData}}}
\examples{
# suite de l'exemple de \code{CECclassif()}
Zmiss <- Z[1:5, ]
Zmiss$X2 <- NA_real_
preds <- CECpredict(Zmiss, fit$params, idColToPred = "X2")
preds
}
