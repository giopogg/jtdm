% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/joint_trait_prob.R
\name{joint_trait_prob}
\alias{joint_trait_prob}
\title{Computes joint probabilities.}
\usage{
joint_trait_prob(
  m,
  indexTrait,
  bounds,
  Xnew = NULL,
  FullPost = T,
  mcmc.samples = NULL,
  parallel = FALSE
)
}
\arguments{
\item{m}{a model fitted with \code{jtdm_fit}}

\item{indexTrait}{A vector of the names (as specified in the column names of Y) of the two (or more!) traits we want to compute the joint probabilities of.}

\item{bounds}{The parameter to specify a region in the community-trait space where the function computes the joint probabilities of traits. It is a list of the length of "indexTrait", each element of the list is a vector of length two. The vector represents the inferior and superior bounds of the region for the specified trait. For example, if we consider two traits, bounds=list(c(10,Inf),c(10,Inf)) corresponds to the region in the community-trait space where both traits both take values greater than 10.}

\item{Xnew}{Optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.}

\item{FullPost}{If FullPost = TRUE, the function returns samples from the predictive distribution of joint  probabilities. If FullPost= FALSE, joint probabilities are computed only using the posterior mean of the parameters. FullPost cannot be equal to "mean" here.}

\item{mcmc.samples}{Optional, default to NULL, only works when FullPost=FALSE. Defines the number of MCMC samples to compute the posterior distribution of joint probabilities. Needs to be between 1 and m$model$sample x length(m$model$mcmc)}

\item{parallel}{Optional, only works when FullPost = TRUE. When TRUE, the function uses mclapply to parallelise the calculation of the posterior distribution joint probabilities.}
}
\value{

}
\description{
Computes the joint probability of CWM traits in regions in the community-trait space specified by bounds and in sites specified in Xnew.
}
\details{
This function is time consuming when \code{FullPost=T}. Consider setting \code{parallel=T} and/or to set \code{mcmc.samples} to a value smaller than the length of the MCMC chains.
}
\examples{
data(Y)  
data(X)  
# Short MCMC to obtain a fast example: results are unreliable !
m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10,
        burnin = 100,
        sample = 100)  
# Compute probability of SLA and LNC to be joint-high at sites in the studies
joint = joint_trait_prob(m,indexTrait=c("SLA","LNC"),
                         bounds=list(c(mean(Y[,"SLA"]),Inf),c(mean(Y[,"SLA"]),Inf)))
}
