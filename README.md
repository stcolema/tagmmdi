
<!-- README.md is generated from README.Rmd. Please edit that file -->
TAGMMDI
============

[![Travis build status](https://travis-ci.org/stcolema/BayesicGibbs.svg?branch=master)](https://travis-ci.org/stcolema/BayesicGibbs) [![AppVeyor build status](https://ci.appveyor.com/api/projects/status/github/stcolema/BayesicGibbs?branch=master&svg=true)](https://ci.appveyor.com/project/stcolema/BayesicGibbs) [![Coverage status](https://codecov.io/gh/stcolema/BayesicGibbs/branch/master/graph/badge.svg)](https://codecov.io/github/stcolema/BayesicGibbs?branch=master)

The goal of tagmmdi is to extend Crook's T-Augmented Gaussian Mixture (TAGM) model to incorporate Multiple Dataset Integration (MDI).

Installation
------------

# You can install the released version of tagmmdi from [CRAN](https://CRAN.R-project.org) with:
You can install the package from GITHUB using
``` r
devtools::install_github("tagmmdi")
```

Example
-------

This is a basic example of applying MDI clustering to data from pRolocData:

``` r


```

What is special about using `README.Rmd` instead of just `README.md`? You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don't forget to commit and push the resulting figure files, so they display on GitHub!
