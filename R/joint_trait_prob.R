#' Computes joint probabilities.
#'
#' Computes the joint probability of CWM traits in regions in the community-trait space specified by bounds and in sites specified in Xnew.
#' @param m a model fitted with \code{jtdm_fit}
#' @param indexTrait  A vector of the names (as specified in the column names of Y) of the two (or more!) traits we want to compute the joint probabilities of.
#' @param bounds The parameter to specify a region in the community-trait space where the function computes the joint probabilities of traits. It is a list of the length of "indexTrait", each element of the list is a vector of length two. The vector represents the inferior and superior bounds of the region for the specified trait. For example, if we consider two traits, bounds=list(c(10,Inf),c(10,Inf)) corresponds to the region in the community-trait space where both traits both take values greater than 10.
#' @param Xnew Optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used.
#' @param FullPost If FullPost = TRUE, the function returns samples from the predictive distribution of joint  probabilities. If FullPost= FALSE, joint probabilities are computed only using the posterior mean of the parameters. FullPost cannot be equal to "mean" here.
#' @param mcmc.samples Optional, default to NULL, only works when FullPost=FALSE. Defines the number of MCMC samples to compute the posterior distribution of joint probabilities. Needs to be between 1 and m$model$sample x length(m$model$mcmc)
#' @param parallel Optional, only works when FullPost = TRUE. When TRUE, the function uses mclapply to parallelise the calculation of the posterior distribution joint probabilities.
#' @export
#' @return A list containing:
#'    \item{PROBsamples}{Samples from the posterior distribution of the joint probability.NULL if FullPost=FALSE. }
#'    
#'    \item{PROBmean}{Posterior mean of the joint probability.}
#'    
#'    \item{PROBq975,PROBq025}{97.5\% and 0.25\% posterior quantiles of the joint probability. NULL if FullPost=FALSE. }
#' }
#' @details This function is time consuming when \code{FullPost=T}. Consider setting \code{parallel=T} and/or to set \code{mcmc.samples} to a value smaller than the length of the MCMC chains.
#' @examples
#' data(Y)  
#' data(X)  
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10,
#'         burnin = 100,
#'         sample = 100)  
#' # Compute probability of SLA and LNC to be joint-high at sites in the studies
#' joint = joint_trait_prob(m,indexTrait=c("SLA","LNC"),
#'                          bounds=list(c(mean(Y[,"SLA"]),Inf),c(mean(Y[,"SLA"]),Inf)))
joint_trait_prob = function(m, indexTrait, bounds, Xnew=NULL, FullPost=T, mcmc.samples=NULL, parallel=FALSE){

  ntot = m$model$sample*length(m$model$mcmc) #samples * n.chains
  if(!is.null(mcmc.samples)){if(mcmc.samples>ntot){stop("You need to provide a number of mcmc samples lower than the length of the chain given by m$model$sample*length(m$model$mcmc)")}}
  if(is.null(Xnew)) Xnew=m$X_raw
  if(is.null(dim(Xnew))) Xnew=t(as.matrix(Xnew))
  if(length(indexTrait) != length(bounds)){stop("index and bounds have different lengths!!")}
  if(ncol(Xnew)!=ncol(m$X_raw)) stop("provide Xnew with the same number of X")
  data=list(Y=m$Y, X=m$X, K=ncol(m$X), J=ncol(m$Y), n=nrow(m$Y), df= ncol(m$Y), I=diag(ncol(m$Y)),  X_raw = m$X_raw, poly_coef= m$poly_coef, poly.var = m$poly.var, poly.degree = m$poly.degree)

  if(!identical(colnames(Xnew),colnames(m$X_raw)))  stop("Provide Xnew with the same colnames of the design matrix used to fit the model")

  index = sapply(indexTrait,function(x){which(colnames(m$Y) %in% x )})

  Xnew=model.frame(m$mt,as.data.frame(Xnew))
  Xnew=model.matrix(m$mt,Xnew)

  Xnew=apply(Xnew,MARGIN=c(1,2),FUN=as.numeric)

  lower = sapply(bounds, "[[", 1)
  upper = sapply(bounds, "[[", 2)

  if(FullPost){

    B=getB(m)$Bsamples
    Sigma=get_sigma(m)$Ssamples
    ntot = m$model$sample*length(m$model$mcmc) #samples * n.chains

    PROB = array(dim=c(nrow(Xnew),ntot))

    if(!is.null(mcmc.samples)){
      if(parallel){
        PROB=parallel::mclapply(sample(1:ntot,mcmc.samples,replace=F),FUN=function(i){
          apply(Xnew %*% t(B[index,,i]),MARGIN=1,FUN=function(x)
            pmvnorm(lower, upper, mean=x,  sigma=Sigma[index,index,i]) )},  mc.cores = parallel::detectCores(), mc.allow.recursive = TRUE)
        PROB = matrix(unlist(PROB),ncol=mcmc.samples)
    }else{
      PROB=sapply(sample(1:ntot,mcmc.samples,replace=F),FUN=function(i){
        apply(Xnew %*% t(B[index,,i]),MARGIN=1,FUN=function(x)
          pmvnorm(lower, upper, mean=x,  sigma=Sigma[index,index,i]) )})
    }
    }else{
      if(parallel){
        PROB=parallel::mclapply(1:ntot,FUN=function(i){
          apply(Xnew %*% t(B[index,,i]),MARGIN=1,FUN=function(x)
            pmvnorm(lower, upper, mean=x,  sigma=Sigma[index,index,i]) )},  mc.cores = parallel::detectCores(), mc.allow.recursive = TRUE)
        PROB = matrix(unlist(PROB),ncol=mcmc.samples)
      }else{
     PROB=sapply(1:ntot,FUN=function(i){
       apply(Xnew %*% t(B[index,,i]),MARGIN=1,FUN=function(x)
         pmvnorm(lower, upper, mean=x,  sigma=Sigma[index,index,i]) )})
    }
    }
    if(is.null(dim(PROB))){PROB=t(as.matrix(PROB))}
    PROB_hat = apply(PROB,mean,MARGIN = 1)
    PROB_975 = apply(PROB, quantile, MARGIN=1,0.975)
    PROB_025 = apply(PROB, quantile, MARGIN=1,0.025)

    if(!is.null(rownames(Xnew))) names( PROB_hat)=names(PROB_975)=names(PROB_025)=rownames(Xnew)


  }else{

    B=getB(m)$Bmean
    Sigma=get_sigma(m)$Smean
    mu = Xnew %*% t(B[index,])
    varcov = Sigma[index,index]

    PROB_hat = apply(mu,MARGIN=1,FUN=function(x)pmvnorm(lower, upper, mean=x,  sigma=varcov))

    PROB = PROB_975 = PROB_025 = NULL
    if(!is.null(rownames(Xnew))) names(PROB_hat) = rownames(Xnew)
  }
  list(PROBsamples = PROB, PROBmean=PROB_hat, PROBq025=PROB_025,PROBq975=PROB_975)

}
