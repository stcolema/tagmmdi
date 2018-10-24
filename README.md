
<!-- README.md is generated from README.Rmd. Please edit that file -->
tagmmdi
============

[![Travis build status](https://travis-ci.org/stcolema/BayesicGibbs.svg?branch=master)](https://travis-ci.org/stcolema/BayesicGibbs) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/stcolema/BayesicGibbs?branch=master&svg=true)](https://ci.appveyor.com/project/stcolema/BayesicGibbs) [![Coverage status](https://codecov.io/gh/stcolema/BayesicGibbs/branch/master/graph/badge.svg)](https://codecov.io/github/stcolema/BayesicGibbs?branch=master)

The goal of tagmmdi is to extend Crook's T-Augmented Gaussian Mixture (TAGM) model to incorporate Multiple Dataset Integration (MDI).

Installation
------------

You can install the package from GITHUB using
``` r
devtools::install_github("tagmmdi")
```

Example
-------

This is a basic example of applying MDI clustering to data from pRolocData:

``` r
data(tan2009r1)
data(tan2009r1goCC)
cat_data <- as.matrix(exprs(tan2009r1goCC))

mcmc_obj <- mcmc_out(tan2009r1,
                     cat_data = cat_data,
                     num_clusters_cat = 30,
                     num_iter = 5000,
                     burn = 1000,
                     thinning = 25,
                     outlier = T)

```
