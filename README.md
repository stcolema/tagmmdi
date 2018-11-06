
<!-- README.md is generated from README.Rmd. Please edit that file -->
tagmmdi
============

[![Travis build status](https://travis-ci.org/stcolema/tagmmdi.svg?branch=master)](https://travis-ci.org/stcolema/tagmmdi) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/stcolema/tagmmdi?branch=master&svg=true)](https://ci.appveyor.com/project/stcolema/tagmmdi) [![Coverage status](https://codecov.io/gh/stcolema/tagmmdi/branch/master/graph/badge.svg)](https://codecov.io/github/stcolema/tagmmdi?branch=master)

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
data(tan2009r1)
data(tan2009r1goCC)
cat_data <- as.matrix(exprs(tan2009r1goCC))

mcmc_obj <- mcmc_out(tan2009r1,
                     data_2 = cat_data,
                     n_clust_2 = 30,
                     num_iter = 5000,
                     burn = 1000,
                     thinning = 25,
                     outlier_1 = T)

```
