COMBO
==================================================
[![Downloads](http://cranlogs.r-pkg.org/badges/rddapp)](https://CRAN.R-project.org/package=rddapp) [![Cran Build](https://www.r-pkg.org/badges/version/rddapp)](https://CRAN.R-project.org/package=rddapp)
[![R-CMD-check](https://github.com/felixthoemmes/rddapp/actions/workflows/r.yml/badge.svg)](https://github.com/felixthoemmes/rddapp/actions/workflows/r.yml)

Overview
--------------------------------------------------

**COMBO** provides a set of functions for the analysis of regression models with binary outcome misclassification. 

The two main parts are:

- Classification probability calculations
- Parameter estimation 


Classification probability calculations
--------------------------------------------------
The package allows users to compute the probability of the latent true outcome and the conditional probability of observing an outcome given the latent true outcome, based on parameters estimated from the `COMBO_EM` and `COMBO_MCMC` functions.


Parameter estimation 
--------------------------------------------------
Jointly estimate parameters from the true outcome and observation mechanisms, respectively, in a binary outcome misclassification model using the EM algorithm or MCMC. 

Installation
--------------------------------------------------

``` r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("kimhochstedler/COMBO")
```
