
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)

# jtdm - Joint trait distribution modeling

## About the method

The package jtdm implements the method described in Poggiato et al., in
prep. Joint models and predictions of community traits. The code for
producing the results in this paper is available under the subfolder
publications in this repo.

## Installing the R package

### R-package

The package implements jtdm using the Markov Chain Monte Carlo Bayesian
modeling software JAGS via the R package runjags. Therefore, it requires
the installation of JAGS. Its installation is easy and depends on your
operating system:

##### Ubutntu

sudo apt-get install jags

##### Windows

<https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/>

##### Mac

<https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/>

Once JAGS has been installed, the following code should run:

``` r
library(devtools)
install_github("giopogg/jtdm")
```

The package implements jtdm using the Markov Chain Monte Carlo Bayesian
modeling software JAGS via the R package runjags. Therefore, it requires
the installation of both JAGS and runjags.

Fitting a JTDM the case study dataset of Poggiato et al. In prep.

``` r
library(jtdm)
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
    ## SLA       8.008574 0.010900906 0.4912318 10.621769
    ## LNC      20.383691 0.001789419 0.1440694  2.889193
    ## Height   10.369534 0.018586089 0.1465588  2.981464

``` r
get_sigma(m)$Smean
```

    ##              SLA       LNC    Height
    ## SLA     77.05206 17.760589 -12.61709
    ## LNC     17.76059  9.261374   1.04773
    ## Height -12.61709  1.047730  94.63345

Trait-environment relationships

``` r
 partial_response(m,indexGradient="GDD",indexTrait="SLA",FixX=list(GDD=NULL,FDD=NULL,forest=1))$p
```

Partial response curve of the most suitable community-level strategy and
envelop of possible community-level strategies of SLA and LNC along the
GDD gradient.

``` r
ellipse_plot(m,indexTrait = c("SLA","LNC"),indexGradient="GDD")
```

![](README_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

Computes joint probabilities of both SLA and LNC to be greater than 20
in a high altitude site. This measures the relative suitability of
communities where both SLA and LNC are higher than 20 in a high altitude
site

``` r
joint_trait_prob(m,indexTrait=c("SLA","LNC"), Xnew=X["VCHA_2940",], bounds=list(c(20,Inf),c(20,Inf)))$PROBmean
```

    ##         1 
    ## 0.1057429

Unsurprisingly, the probability is low. Then, we compute this
probability along the GDD gradient

``` r
joint=joint_trait_prob_gradient(m,indexTrait=c("SLA","LNC"), indexGradient="GDD", bounds=list(c(mean(Y[,"SLA"]),Inf),c(mean(Y[,"SLA"]),Inf)))
```

And plot it

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- --> As climatic
conditions become more favourable (i.e. GDD increases), the probability
of having high values of both traits increases.
