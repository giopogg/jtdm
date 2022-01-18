
<!-- README.md is generated from README.Rmd. Please edit that file -->

``` r
knitr::opts_chunk$set(fig.path='Figs/')
```

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# jtdm - Joint trait distribution modeling

## About the method

The method is described in Poggiato et al., in prep. Joint models and
predictions of community traits. The code for producing the results in
this paper is available under the subfolder publications in this repo.

The method itself is implemented in this repo.

## Installing the R package

### R-package

Install the package via

``` r
install_github("giopogg/jtdm")
```

The package implements jtdm using the Markov Chain Monte Carlo Bayesian
modeling software JAGS via the R package runjags. Therefore, it requires
the installation of both JAGS and runjags.

Once the dependencies are installed, the following code should run:

Fits a JTDM

``` r
library(jtdm)
```

    ## Warning: remplacement de l'importation précédente 'arm::traceplot' par
    ## 'coda::traceplot' lors du chargement de 'jtdm'

``` r
library(ggplot2)
set.seed(1712)
data(Y)
data(X)
# Short MCMC to obtain a fast example: results are unreliable !
m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10, burnin = 100, sample = 100)
```

    ## module dic loaded

    ## Compiling rjags model...
    ## Calling the simulation using the rjags method...
    ## Note: the model did not require adaptation
    ## Burning in the model for 100 iterations...
    ## Running the model for 100 iterations...
    ## Extending 100 iterations for pD/DIC estimates...
    ## Simulation complete
    ## Calculating summary statistics...
    ## Calculating the Gelman-Rubin statistic for 21 variables....
    ## Note: Unable to calculate the multivariate psrf
    ## Finished running the simulation

``` r
# Inferred parameters
getB(m)$Bmean
```

    ##        (Intercept)         GDD       FDD    forest
    ## SLA       10.50979 0.008854437 0.4896142 11.769861
    ## LNC       21.39144 0.001005003 0.1472568  3.366346
    ## Height    12.13099 0.016846515 0.1165499  4.043582

``` r
get_sigma(m)$Smean
```

    ##              SLA        LNC      Height
    ## SLA     75.46582 17.1209341 -13.7401827
    ## LNC     17.12093  9.0821226   0.9269098
    ## Height -13.74018  0.9269098  95.7111518

Computes joint probabilities of both SLA and LNC to be greater than 10
in the sites of the study.

``` r
joint = joint_trait_prob(m,indexTrait=c("SLA","LNC"), bounds=list(c(mean(Y[,"SLA"]),Inf),c(mean(Y[,"SLA"]),Inf)))
```
