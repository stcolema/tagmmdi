
<!-- README.md is generated from README.Rmd. Please edit that file -->
tagmmdi
============

The goal of tagmmdi is to extend Crook's T-Augmented Gaussian Mixture (TAGM) model to incorporate Multiple Dataset Integration (MDI).

Installation
------------

You can install the package from GITHUB using
``` r
devtools::install_github("tagmmdi")
```

Example
-------

This is a basic example of applying MDI clustering to data from pRolocdata:

``` r
# Import data from pRolocdata
data(tan2009r1)
data(tan2009r1goCC)

# Run MDI clustering on Gaussian and Categorical datasets
mcmc_obj <- mcmc_out(tan2009r1,
                     data_2 = as.matrix(exprs(tan2009r1goCC)),
                     n_clust_2 = 30,
                     num_iter = 5000,
                     burn = 1000,
                     thinning = 25,
                     outlier_1 = T)

```
