# CEClust

**Composite Entropy Clustering for mixed data (R package)**

[![R-CMD-check](https://github.com/dumont-thierry/CEClust/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/dumont-thierry/CEClust/actions)

## Overview

**CEClust** performs unsupervised clustering on **mixed-type data** (numeric + categorical) using the **Composite Entropy Criterion (CEC)**.  
It learns clusters, estimates parameters, and can predict/complete missing columns conditional on the inferred clusters.

Core functions:

- `CECclassif()` – fit the clustering on a dataset
- `CECclassifNewData()` – assign new observations to the learned clusters
- `CECpredict()` – predict missing columns given the cluster structure

## Installation

Install the development version from GitHub:

```r
# Install remotes if needed
install.packages("remotes")

# Install CEClust from GitHub
remotes::install_github("dumont-thierry/CEClust")
```

Load the package:

```r
library(CEClust)
```

## Quick example

```r
set.seed(123)
n <- 100

# Two Gaussian clusters on 2 variables
x1 <- rnorm(n, 0, 1); x2 <- rnorm(n, 0, 1)
x3 <- rnorm(n, 3, 1); x4 <- rnorm(n, 3, 1)
X  <- rbind(cbind(x1, x2), cbind(x3, x4))

# A categorical variable correlated with clusters
lab <- factor(rep(c("A","B"), each = n))
Z   <- data.frame(X1 = X[,1], X2 = X[,2], Cat = lab)

# Fit
fit <- CECclassif(Z, Nshots = 50, Nloop = 200)

# Cluster assignments
table(fit$clusters, lab)

# Assign new observations (here: only numeric columns provided)
Znew <- Z[, c("X1","X2")]
pred_clusters <- CECclassifNewData(Znew, fit$params, idColToPred = "Cat")

# Predict a missing numeric column given the cluster
Znew2 <- Z[, c("X2","Cat")]
pred_x1 <- CECpredict(Znew2, fit$params, idColToPred = "X1")
```

## Documentation

After installation, run in R:

```r
?CECclassif
?CECclassifNewData
?CECpredict
```

These help pages include arguments, return values, references to sibling functions, and runnable examples.

## Dependencies

- Base R (`stats`)
- `mvtnorm`

## License

MIT License. See `LICENSE`.
