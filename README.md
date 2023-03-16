# COMBO ![my workflow](https://github.com/kimhochstedler/COMBO/actions/workflows/r.yml/badge.svg)

![ ](https://github.com/kimhochstedler/COMBO/blob/main/small_logo.png?raw=true)

**COMBO:** **CO**rrecting **M**isclassified **B**inary **O**utcomes

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
Jointly estimate parameters from the true outcome and observation mechanisms, respectively, in a binary outcome misclassification model using the EM algorithm or MCMC. Parameters from the true outcome, first-stage observation, and second-stage observation mechanisms in a two-stage binary outcome misclassification model can also be estimated using the EM algorithm and MCMC.

Installation
--------------------------------------------------

``` r
# Install the development version from GitHub:
# install.packages("devtools")
devtools::install_github("kimhochstedler/COMBO")
```
