
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
the installation of JAGS.

##### Ubutntu

sudo apt-get install jags

##### Windows

<https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Windows/>

##### Mac

<https://sourceforge.net/projects/mcmc-jags/files/JAGS/4.x/Mac%20OS%20X/>

Once JAGS have been installed, the following code should run:

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
    ## SLA        8.78095 0.010872252 0.5490040 11.059551
    ## LNC       20.65485 0.001937326 0.1814652  3.021635
    ## Height    12.36151 0.016775119 0.1248564  3.824376

``` r
get_sigma(m)$Smean
```

    ##              SLA        LNC      Height
    ## SLA     75.81373 17.3085241 -14.5933684
    ## LNC     17.30852  9.1407323   0.6651481
    ## Height -14.59337  0.6651481  95.4718023

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
in a high altitude site.

``` r
joint_trait_prob(m,indexTrait=c("SLA","LNC"), Xnew=X["VCHA_2940",], bounds=list(c(20,Inf),c(20,Inf)))$PROBmean
```

    ##         1 
    ## 0.1077798

Then, we compute this probability along the GDD gradient

``` r
joint=joint_trait_prob_gradient(m,indexTrait=c("SLA","LNC"), indexGradient="GDD", bounds=list(c(mean(Y[,"SLA"]),Inf),c(mean(Y[,"SLA"]),Inf)))
```

And plot it

![](README_files/figure-gfm/unnamed-chunk-6-1.png)<!-- -->
