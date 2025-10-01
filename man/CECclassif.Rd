\name{CECclassif}
\alias{CECclassif}
\title{Clustering par Entropie Composite sur données mixtes}
\description{
Apprend un modèle de mélange (gaussien pour les variables numériques et multinomial pour les variables catégorielles) et renvoie l’assignation de cluster ainsi que les paramètres estimés.
}
\usage{
CECclassif(Z, lambda = 1, C = 1, r0 = NULL, Nshots = 100, Nloop = 1000,
           familyType = "gaussAndDiscreteVector", sizeMaxOutlier = 0,
           autoRegroupOutliers = FALSE, displayRemainingTime = FALSE, focus = NULL)
}
\arguments{
  \item{Z}{\code{data.frame} ou \code{matrix} d'observations (colonnes = variables).}
  \item{lambda}{\code{numeric} (>0). Paramètre d’entropie composite.}
  \item{C}{\code{numeric} (>0). Constante de borne.}
  \item{r0}{\code{integer} ou \code{NULL}. Nombre de classes initial si fourni.}
  \item{Nshots}{\code{integer}. Nombre de redémarrages aléatoires (multi-start).}
  \item{Nloop}{\code{integer}. Max itérations internes par redémarrage.}
  \item{familyType}{\code{character}. Famille cible (\code{"gaussAndDiscreteVector"}, etc.).}
  \item{sizeMaxOutlier}{\code{integer}. Seuil (taille min) d’un groupe considéré outlier à regrouper (si > 0).}
  \item{autoRegroupOutliers}{\code{logical}. Regroupement auto des petits groupes.}
  \item{displayRemainingTime}{\code{logical}. Affiche le temps estimé.}
  \item{focus}{\code{integer} ou \code{character} ou \code{NULL}. Si une colonne factor est fournie, le clustering est appris séparément par niveau puis agrégé.}
}
\value{
Une \code{list} avec au minimum :
\itemize{
  \item \code{phi} : vecteur d'appartenance (one-hot) de longueur \code{n*r}.
  \item \code{params} : liste des paramètres estimés (nu, m, Sigma, etc.).
  \item \code{Hphi} : valeur de l’objectif (entropie composite).
  \item \code{clusters} : vecteur des clusters prédits (longueur \code{n}).
  \item \code{REO} : nombre final de classes retenues.
}}
\seealso{\code{\link{CECclassifNewData}}, \code{\link{CECpredict}}}
\examples{
set.seed(123)
n <- 100
# Deux clusters gaussiens sur 2 variables
x1 <- rnorm(n, 0, 1); x2 <- rnorm(n, 0, 1)
x3 <- rnorm(n, 3, 1); x4 <- rnorm(n, 3, 1)
X  <- rbind(cbind(x1, x2), cbind(x3, x4))
# Une variable catégorielle corrélée au cluster
lab <- factor(rep(c("A","B"), each = n))
Z   <- data.frame(X1 = X[,1], X2 = X[,2], Cat = lab)

fit <- CECclassif(Z, Nshots = 100)
table(fit$clusters, lab)
}
