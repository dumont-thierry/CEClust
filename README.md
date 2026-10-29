# CEClust : grille adaptative et bandes

Version **0.4.0.9000**, assemblée le 29 octobre 2026.

## Installation depuis GitHub

R >= 4.3 et un compilateur C++17 sont nécessaires (Rtools sous Windows).
La construction des guides utilise Pandoc, fourni notamment avec RStudio.

```r
install.packages(c("remotes", "knitr", "rmarkdown"))
remotes::install_github("dumont-thierry/CEClust", build_vignettes = TRUE)
```

## Choisir la méthode

```r
library(CEClust)
fit <- CECdiagram(X, method = "adaptive",
  L = 20, B = 20, workers = 4, seed = 1, directory = "resultats_grille")
fit_bandes <- CECdiagram(X, method = "bands",
  L = 20, B = 20, workers = 4, seed = 1, directory = "resultats_bandes")
```

La grille adaptative est le choix par défaut. Les bandes restent expérimentales.
`penalty = "bic"` active une pénalité d'ordre dans chaque moteur.
Les banques, l'exploration et les règles de validation propres à chaque méthode
restent distinctes. Le choix de méthode ne garantit pas des partitions identiques.

## Grille adaptative

Moteur CEClust 0.3.0.9010 : CEM, propagation et subdivisions conservées.
La nouvelle entrée calcule Cmin comme le pic de la gaussienne ajustée au groupe
unique, avec covariance empirique divisée par n. Les dix centres initiaux valent
`2^(1:10) * Cmin`. Pas d'extension à gauche par défaut. Les extensions à droite
et en lambda sont conservées. Les bords des cellules restent à mi-pas des centres.
Cmin est une référence d'exploration, pas une garantie universelle de saturation.

```r
CECinitialGrid(X)
ctl <- CECdiagramControl("adaptive",
  optimizer_control = list(Rmax = 20, Ninit = 5),
  adaptive_control = list(max_cells = 12000, extend_left = FALSE),
  rsi_threshold = .95)
fit <- CECdiagram(X, "adaptive", control = ctl, directory = "grille_2")
CECsummariseGrid(fit)
CECplotGrid(fit)
```

Une covariance singulière produit une erreur explicite : supprimer les variables
redondantes ou fournir C_grid explicitement, sans régularisation cachée.
Le recalcul conditionnel du RSI reste désactivé et l'accord entre champions
requis par défaut, comme en 9010. Les options historiques restent accessibles.
La variante par banques EM s'obtient avec `control = list(optimization = "em_bank")`.
Les anciens appels CECfitBoundGrid gardent leur ancienne grille par défaut et
leurs quatre extensions. Utiliser CECdiagram pour la nouvelle initialisation.

## Bandes

Moteur CEClustBands 0.1.0.9001, descentes et diagnostics conservés.
Budget de l'entrée commune : ordres pilote 1:20, 5 initialisations, 2 départs
supplémentaires, 10 ancres, ordre enrichi maximal 40, 11/41 colonnes C,
2 tours de raffinement, 2 tours de projection, 200 requêtes.
L=B=20, stabilité > .85, RSI > .95, saturation de rétention nulle.
workers vaut 1 par défaut et doit être choisi explicitement.

```r
ctl <- CECdiagramControl("bands", max_anchors = 10, extra_starts = 2)
fit_bandes <- CECdiagram(X, "bands", control = ctl, directory = "bandes_2")
bands_partitions(fit_bandes)
bands_plot(fit_bandes)
bands_validate(fit_bandes)
# Facultatif, séparé, jamais automatique :
# audit <- bands_cleanup(fit_bandes, directory = "menage", workers = 4)
```

Le ménage élimine seulement : score meilleur ET (entropie des poids plus grande
OU saturation), après subdivision et CEM. Aucun remplacement ni nouveau RSI.
cec_bands et bands_control gardent leurs défauts historiques, plus coûteux.
Ne pas les confondre avec le profil de CECdiagramControl("bands").

## Archives et reprise

Un dossier par essai. Conserver le dossier avec le résultat : les bandes y
référencent leurs allocations. Reprise avec resume=TRUE, mêmes données et
paramètres. Les checkpoints des packages séparés ont un autre runtime et ne
sont pas réutilisés. Les anciennes archives ne sont pas modifiées.

## Sources

- R/diagram.R : entrée commune, Cmin, profils de paramètres.
- R/adaptive_*.R : grille adaptative et géométrie.
- R/bands_*.R : moteur bandes, fichiers nommés selon leur fonction.
- src/bands_moments.cpp : noyau CEM des bandes.

Les sources historiques et le manuscrit soumis restent inchangés.
Ce regroupement ne résout pas les limites d'exploration connues des bandes.

## Versions archivées

- [0.3.0.9002 exacte](https://github.com/dumont-thierry/CEClust/tree/archive/0.3.0.9002),
  commit `61dd3ca85d5f2fe817aa15f8318879c061f750a7`.
- [0.3.0.9008, état de GitHub avant l'unification](https://github.com/dumont-thierry/CEClust/tree/archive/0.3.0.9008),
  commit `2075babd81db884b8906d8c04eeffe537622f832`.

Ces archives conservent les fichiers d'origine. Pour installer la version
0.3.0.9002 exacte :

```r
remotes::install_github("dumont-thierry/CEClust", ref = "archive/0.3.0.9002",
  build_vignettes = TRUE)
```
