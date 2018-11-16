
<!-- README.md is generated from README.Rmd. Please edit that file -->
tagmmdi
============

The goal of tagmmdi is to extend Crook's T-Augmented Gaussian Mixture (TAGM) model to incorporate Multiple Dataset Integration (MDI). Please note that this is still pre-release and probably does not work yet.

Installation
------------

You can install the package from GitHub using
``` r
devtools::install_github("tagmmdi")
```

Initial steps
-------------
The main MDI functions are all designed to run on objects from pRolocdata, so a pre-requisite to running any functions is to load this library and any datasets of interest:

```r
# Load the tagmmdi library to access it
library(tagmmdi)

#Load pRolocdata to access datasets
library(pRolocdata)

# Import data
data(tan2009r1)
data(tan2009r1goCC)

# Optionally use only the informative GO terms to reduce runtime (the 20 here is arbitrary)
useful_GO <- tan2009r1goCC[, colSums(exprs(tan2009r1goCC)) > 20]
```

MDI clustering
--------------

This is a basic example of applying MDI clustering to data from pRolocdata:

``` r
# Run MDI clustering on Gaussian and Categorical datasets
mcmc_obj <- mcmc_out(tan2009r1,
  data_2 = as.matrix(exprs(tan2009r1goCC)),
  n_clust_2 = 30,
  num_iter = 5000,
  burn = 1000,
  thinning = 25,
  outlier_1 = T
)
```

Cross validation
----------------

This is an example of running cross-validation to measure the quadratic loss score of the method:

```r
# Run cross-validation
cv_obj_mdi <- mdi_cross_validate(tan2009r1,
  MS_cat_object = useful_GO,
  times = 10,
  test_size = 0.2,
  num_iter = 5000,
  burn = 1000,
  thinning = 25
)
```
