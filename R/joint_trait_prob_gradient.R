#' Computes partial response curves of joint probabilities
#'
#' Computes the partial responses curves of joint probability of CWM traits as a function of a focal variable. The regions in which joint probabilities are computed are specified by bounds. In order to build the response curve, the function builds a dataframe where the focal variable varies along a gradient and the other (non-focal) variables are fixed to their mean (but see FixX parameter for fixing non-focal variables to user-defined values). Then, uses joint_trait_prob to compute the joint probability in these dataset.
#' @param m A model fitted with \code{jtdm_fit}
#' @param indexTrait  A vector of the names (as specified in the column names of Y) of the two (or more!) traits we want to compute the joint probabilities of.
#' @param indexGradient The name (as specified in the column names of X) of the focal variable.
#' @param bounds The parameter to specify a region in the community-trait space where the function computes the joint probabilities of traits. It is a list of the length of "indexTrait", each element of the list is a vector of length two. The vector represents the inferior and superior bounds of the region for the specified trait. For example, if we consider two traits, bounds=list(c(10,Inf),c(10,Inf)) corresponds to the region in the community-trait space where both traits both take values greater than 10.
#' @param grid.length The number of points along the gradient of the focal variable. Default to 200.
#' @param XFocal Optional. A gradient of the focal variable provided by the user. If provided, the function will used this gradient instead of building a regular one. Default to NULL.
#' @param FixX Optional. A parameter to specify the value to which non-focal variables are fixed. This can be useful for example if we have some categorical variables (e.g. forest vs meadows) and we want to obtain the partial response curve for a given value of the variable. It has to be a list of the length and names of the columns of X. For example, if the columns of X are "MAT","MAP","Habitat" and we want to fix "Habitat" to 1, then FixX=list(MAT=NULL,MAP=NULL,Habitat=1.). Default to NULL.
#' @param FullPost If FullPost = TRUE, the function returns samples from the predictive distribution of joint  probabilities, thus allowing the computation of credible intervals. If FullPost= FALSE, joint probabilities are computed only using the posterior mean of the parameters. FullPost cannot be equal to "mean" here.
#' @param samples Optional, default to NULL, only works when FullPost=FALSE. Defines the number of samples to compute the posterior distribution of joint probabilities. Needs to be between 1 the total number of samples drawn from the posterior distribution.
#' @param parallel Optional, only works when FullPost = TRUE. When TRUE, the function uses mclapply to parallelise the calculation of the posterior distribution joint probabilities.
#' @details This function is time consuming when \code{FullPost = TRUE}. Consider setting \code{parallel = TRUE} and/or to set \code{samples} to a value smaller than the total number of posterior samples.
#' @export
#' @return A list containing:
#'    \item{GradProbssamples}{Sample from the posterior distribution of the joint probability along the gradient. It is a vector whose length is the number of posterior samples. NULL if FullPost=FALSE. }

#'    \item{GradProbsmean}{Posterior mean of the joint probability along the gradient. }
#'    
#'    \item{GradProbsq975,GradProbsq025}{97.5\% and 0.25\% posterior quantiles of the joint probability along the gradient. NULL if FullPost=FALSE. }
#'    
#'    \item{gradient}{The gradient of the focal variable built by the function.}
#' 
#' @examples
#' data(Y)  
#' data(X)  
#' # We sample only few samples from the posterior in order to reduce 
#' # the computational time of the examples.
#' # Increase the number of samples to obtain robust results
#' m = jtdm_fit(Y = Y, X = X, formula = as.formula("~GDD+FDD+forest"),  sample = 10)  
#' # Compute probability of SLA and LNC to be joint-high at sites in the studies
#'
#' # Compute the joint probability of SLA and LNC 
#' # to be joint-high along the GDD gradient
#' joint = joint_trait_prob_gradient(m,indexTrait = c("SLA","LNC"), 
#'                                   indexGradient = "GDD",
#'                                   bounds = list(c(mean(Y[,"SLA"]),Inf),c(mean(Y[,"SLA"]),Inf)),
#'                                   FullPost = TRUE)
#'                                   
#' # Compute the joint probability of SLA and LNC to be joint-high along the
#' # GDD gradient when forest = 1 (i.e. in forests) 
#' joint = joint_trait_prob_gradient(m, indexTrait = c("SLA","LNC"),
#'                                   indexGradient = "GDD",
#'                                   bounds = list(c(mean(Y[,"SLA"]),Inf), c(mean(Y[,"SLA"]),Inf)),
#'                                   FixX = list(GDD = NULL, FDD = NULL, forest = 1),
#'                                   FullPost = TRUE)
#'
#' @importFrom stats quantile 
#' @importFrom parallel mclapply detectCores

joint_trait_prob_gradient = function(m, indexTrait, 
                                     indexGradient, bounds, grid.length = 200, XFocal = NULL,
                                     FixX = NULL, FullPost = FALSE, samples = NULL, parallel = FALSE){
  
  indexTrait = sapply(indexTrait,function(x){which(colnames(m$Y) %in% x )})
  indexGradient = which(colnames(m$X_raw) == indexGradient)
  
  if(!inherits(m, "jtdm_fit")) stop("m is not an object of class jtdm_fit")
  
  ntot = dim(m$model$B)[3] #samples * n.chains
  if(is.null(samples)){ samples = ntot }
  if(samples > ntot){stop("You need to provide a number of samples lower than the number of samples of the posterior distribution (parameter `sample` used in `jtdm_fit`")}
  if(FullPost == "mean"){stop("Fullstop cannot be set to mean here.")}
  if(length(indexTrait) != length(bounds)){stop("index and bounds have different lengths!!")}
  if(length(indexGradient) > 1) {stop("indexGradient has to be a numerical value")}
  if(!is.null(FixX)){if(!identical(names(FixX), colnames(m$X_raw))) {stop("Provide FixX as a list with the same names of X")}}
  if(!is.null(FixX[[names(FixX)[indexGradient]]])) {stop("FixX of the focal environmental variable must be NULL")}
  
  data=list(Y = m$Y, X = m$X, K = ncol(m$X), J = ncol(m$Y), n = nrow(m$Y), df = ncol(m$Y),
            I = diag(ncol(m$Y)), X_raw = m$X_raw)
  
  
  #Build XGradient matrix
  if(is.null(XFocal)){
    XGradientFocal = seq(from=  min(data$X_raw[,indexGradient]), to = max(data$X_raw[,indexGradient]), 
                         length.out = grid.length)
  }else{
    XGradientFocal = XFocal
    grid.length = nrow(XFocal)
  }
  
  XGradient = matrix(nrow = grid.length, ncol = ncol(data$X_raw))
  for(i in 1:ncol(data$X_raw)){
    if(i == indexGradient){ XGradient[,i] = XGradientFocal
    }else{
      if(!is.null(FixX[[colnames(data$X_raw)[i]]])){
        XGradient[,i] = rep(FixX[[colnames(data$X_raw)[i]]],grid.length)
      }else{
        XGradient[,i] = rep(mean(data$X_raw[,i]),grid.length)
      }
    }
  }
  colnames(XGradient) = colnames(data$X_raw)
  
  GradProbs = matrix(nrow=grid.length,ncol=ntot)
  # Loop on every row of X
  if(FullPost){
    
    GradProbs = joint_trait_prob(m, indexTrait = colnames(m$Y)[indexTrait],
                                 Xnew = XGradient, bounds = bounds, FullPost = TRUE, 
                                 samples = samples, parallel = parallel)$PROBsamples
    GradProbs_hat = apply(GradProbs, mean, MARGIN=1)
    GradProbs_975 = apply( GradProbs, quantile, MARGIN=1,0.975)
    GradProbs_025 = apply( GradProbs, quantile, MARGIN=1,0.025)
    
  }else{
    
    GradProbs_hat = joint_trait_prob(m, indexTrait = colnames(m$Y)[indexTrait],
                                     Xnew = XGradient, bounds = bounds, FullPost = FALSE)$PROBmean
    GradProbs_975 = GradProbs_025 = NULL
  }
  

  list(GradProbssamples = GradProbs, GradProbsmean=GradProbs_hat,
       GradProbsq025 = GradProbs_025,GradProbsq975 = GradProbs_975, gradient = XGradientFocal)
  
}
