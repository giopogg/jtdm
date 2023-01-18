#' Fitting joint trait distribution models
#'
#' jtdm_fit is used to fit a Joint trait distribution model. Requires the response variable Y (the sites x traits matrix) and the explanatory variables X.
#' @param Y The sites x traits matrix containing community (weighted) means of each trait at each site.
#' @param X The design matrix, i.e. sites x predictor matrix containing the value of each explanatory variable (e.g. the environmental conditions) at each site.
#' @param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under ‘Details’
#' @param adapt,burnin,sample,n.chains,monitor Parameters of the MCMC sampler. See \code{?run.jags} for details
#' @export
#' @details A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See formula for more details of allowed formulae.
#' @return A list containing:
#'    \item{model}{ An object of class 'runjags' containing the fitted model.}
#'    \item{Y}{A numeric vector of standard errors on parameters}
#'    
#'    \item{X_raw}{The design matrix specified as input}
#'    
#'    \item{X}{The design matrix transformed as specified in formula}
#'    
#'    \item{formula}{The formula specified as input}
#' 
#' @examples
#' data(Y)  
#' data(X)  
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10, 
#'         burnin = 100,  
#'         sample = 100)  
#' @importFrom stats model.frame model.matrix rWishart coef
#' @importFrom runjags run.jags
#' @importFrom arm bayesglm
jtdm_fit = function(Y, X, # ! X must not contain the intercept column too !,
                    formula,
                    adapt = 200,
                    burnin = 5000,
                    sample = 5000,
                    n.chains=2,
                    monitor = c('B','Sigma','pd')){
  
  
  X_raw = X
  X=model.frame(formula,as.data.frame(X))
  mt <- attr(X, "terms")
  X <- model.matrix(mt, X)
  
  # data preparation
  data=list(Y=Y, X=X, K=ncol(X), J=ncol(Y), n=nrow(Y), df= ncol(Y), I=diag(ncol(Y)))
  
  # define jags model
  jags.model="model {
      for (i in 1:n) {
      Y[i, 1:J] ~ dmnorm(Mu[i, ], Tau)
      for (j in 1:J) {
      Mu[i, j] <- inprod(B[j, ], X[i, ])
      }
      }
      Tau ~ dwish(I, df)
      Sigma <- inverse(Tau)
      for(j in 1:J){
      for (k in 1:K) {
      B[j,k] ~ dnorm(0, 1.0E-4)
      }
      }
  }"
  
  # set initial values
  inits <- function(data) {
    Y <- as.matrix(data$Y)
    X <- as.matrix(data$X)[, -1]
    Tau <- rWishart(1, data$df, data$I)[, , 1]
    B <- sapply(
      seq_len(data$J),
      function(x) coef(bayesglm(Y[, x] ~ X, family = "gaussian"))
    )
    B <- t(B)
    list(Tau = Tau, B = B)
  }
  
  initsList = rep(list(inits(data)),n.chains)
  # run model
  jags.out <- runjags::run.jags(jags.model,
                                data = data,
                                inits=initsList,
                                #inits=inits(data),
                                adapt = adapt,
                                burnin = burnin,
                                sample = sample,
                                n.chains=n.chains,
                                monitor = monitor,
                                summarise = TRUE)
  
  list(model=jags.out, Y=Y, X=X, X_raw = X_raw, formula=formula,mt=mt)
  
}





#' Predict method for joint trait distribution model
#'
#' Obtains predictions from a fitted joint trait distribution model and optionally computes their R squared and root mean square error (RMSE)
#' @param m a model fitted with \code{jtdm_fit}
#' @param Xnew optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used
#' @param validation  boolean parameter to decide whether we want to compute goodness of fit measures. If true, then Ynew is needed.
#' @param Ynew   Optional. The observed response variables at sites specified in Xnew. It is used to compute goodness of fit metrics when validation= T.
#' @param FullPost   The type of predictions to be obtain. If FullPost = TRUE, the function returns samples from the predictive distribution. If FullPost="mean", the function computes the posterior distribution of the regression term \eqn{BXnew}). If FullPost=F, the function only returns the posterior mean of the regression term (\eqn{BmeanXnew}). 
#' @details To obtain a full assesment of the posterior distribution, the function should be ran with FullPost=TRUE, altough this can be time consuming. FullPost="mean" is used to compute partial response curves, while FullPost=FALSE is used to compute goodness of fit metrics.
#' @export
#' @return A list containing:
#'    \item{Pred}{Sample from the posterior distribution of the posterior predictive distribution. It is an array where the first dimension is the number of sites in Xnew, the second is the number of traits modelled and the third the number of MCMC samples. NULL if FullPost=FALSE.}
#'    
#'    \item{PredMean}{Posterior mean of posterior predictive distribution }
#'    
#'    \item{Predq975,Predq025}{97.5\% and 0.25\% posterior quantiles of the posterior predictive distribution. NULL if FullPost=FALSE. }
#'    
#'    \item{R2}{R squared of predictions (squared Pearson correlation between Ynew and the predictions). NULL if validation=FALSE. }
#'    
#'    \item{RMSE}{Root square mean error between  squared of predictions. NULL if validation=FALSE.}
#' @examples
#' data(Y)  
#' data(X)  
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10,
#'         burnin = 100,
#'         sample = 100)
#' # marginal predictions of traits in the sites of X
#' pred = jtdm_predict(m)
#' @importFrom stats model.frame model.matrix quantile cor 
#' @importFrom coda as.mcmc
#' @importFrom MASS mvrnorm
jtdm_predict = function(m=m, Xnew=NULL, Ynew = NULL, validation = F, FullPost=T){
  
  data=list(Y=m$Y, X=m$X, K=ncol(m$X), J=ncol(m$Y), n=nrow(m$Y), df= ncol(m$Y), I=diag(ncol(m$Y)),  X_raw = m$X_raw)
  
  if(is.null(Xnew)) Xnew=m$X_raw
  if(is.null(dim(Xnew))){Xnew=t(as.matrix(Xnew))}
  if(validation == T & is.null(Ynew)) stop(" if validation = T, you need to provide Ynew!")
  if(ncol(Xnew) != ncol(data$X_raw)) stop("The number of columns of X and Xnew differ")
  if(!is.null(Ynew)){
    if(nrow(Xnew) != nrow(Ynew)) stop("The number of line of Xnew and Ynew differ")
    if(ncol(Ynew) != ncol(data$Y)) stop("The number of columns of Y and Ynew differ")
  }
  if(is.null(FullPost)) stop("please tell whether you want predictions or fitted")
  if(!identical(colnames(Xnew),colnames(m$X_raw))) stop("Provide same column names and same order of the colums in Xnew!")
  
  ###### trasform Xnew with formula
  Xnew_raw = Xnew
  Xnew=model.frame(m$mt,as.data.frame(Xnew))
  Xnew=model.matrix(m$mt,Xnew)
  
  ### Compute predictions
  if(FullPost != FALSE){
    mcmc_param=suppressWarnings(coda::as.mcmc(m$model))
    
    B=getB(m)$Bsamples
    Sigma=get_sigma(m)$Ssamples
    ntot = m$model$sample*length(m$model$mcmc) #samples * n.chains
    
    Predictions = array(dim=c(nrow(Xnew),ncol(data$Y),ntot))
    for(i in 1:ntot){
      if(FullPost==T){
        Predictions[,,i] = t(apply(Xnew, FUN=function(x) {mvrnorm(n=1, mu=B[,,i]%*%x , Sigma=Sigma[,,i])},MARGIN = 1))
      }else{
        Predictions[,,i] = t(apply(Xnew, FUN=function(x) {B[,,i]%*%x},MARGIN = 1))
      }
    }
    
    meanPred = apply(Predictions,mean,MARGIN = c(1,2))
    Pred975 = apply(Predictions, quantile, MARGIN=c(1,2),0.975)
    Pred025 = apply(Predictions, quantile, MARGIN=c(1,2),0.025)
    colnames(meanPred)=colnames(Pred975)=colnames(Pred025)=colnames(data$Y)
    if(!is.null(rownames(Xnew))) rownames(meanPred)=rownames(Pred975)=rownames(Pred025)=rownames(Xnew)
    
  }else{ #If we only want to compoute the mean
    meanPred = as.matrix(Xnew) %*% t(as.matrix(getB(m)$Bmean))
    
    colnames(meanPred)=colnames(data$Y)
    if(!is.null(rownames(Xnew))) rownames(meanPred)=rownames(Xnew)
    
    Pred025 = Pred975 = NULL
    Predictions = NULL
  }
  
  R2_mod=RMSE_mod=NULL
  if(validation){
    #RMSE
    RMSE = function(obs,pred){sqrt(mean((obs - pred)^2))}
    RMSE_mod = vector()
    for(j in 1:ncol(Ynew)) RMSE_mod[j] = RMSE(Ynew[,j],meanPred[,j])
    names(RMSE_mod)=colnames(data$Y)
    #R2
    R2= function(obs,pred){cor(obs, pred) ^ 2}
    R2_mod = vector()
    for(j in 1:ncol(Ynew)) R2_mod[j] = R2(Ynew[,j],meanPred[,j])
    names(R2_mod)=colnames(data$Y)
  }
  
  list(Pred = Predictions, PredMean = meanPred, Predq025 = Pred025, Predq975 = Pred975, R2 = R2_mod, RMSE = RMSE_mod)
  
}





#' Get the inferred regression coefficients
#'
#' Get the samples from the posterior distribution of the regression coefficient matrix B, together with the posterior mean and quantiles. The regression coefficient matrix B is a matrix where the number of rows is defined by the number of traits that are modelled, and the number of columns is the number of columns of the matrix m$X (the number of explanatory variables after transformation via formula)
#' @param m a model fitted with \code{jtdm_fit}
#' @export
#' @return A list containing:
#'    \item{Bsamples}{Sample from the posterior distribution of the regression coefficient matrix. It is an array where the first dimension is the number of traits, the second the number of columns in m$X (the number of variables after transformation via formula) and the third the number of MCMC samples.}
#'    \item{Bmean}{Posterior mean of the regression coefficient matrix.}
#'    \item{Bq975,Bq025}{97.5\% and 0.25\% posterior quantiles of the regression coefficient matrix.}
#' @examples
#' data(Y)  
#' data(X)
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10, 
#'         burnin = 100, 
#'         sample = 100) 
#' # get the inferred regression coefficients
#' B=getB(m)
#' @importFrom stats quantile
#' @importFrom coda as.mcmc

getB=function(m){
  
  data=list(Y=m$Y, X=m$X, K=ncol(m$X), J=ncol(m$Y), n=nrow(m$Y), df= ncol(m$Y), I=diag(ncol(m$Y)))
  
  mcmc_param=suppressWarnings(coda::as.mcmc(m$model))
  ntot = m$model$sample*length(m$model$mcmc) #samples * n.chains
  
  B= array(dim=c(data$J,data$K,ntot))
  for(i in 1:ntot){
    B[,,i]=matrix(mcmc_param[i,grep("B",colnames(mcmc_param))],ncol=data$K)
  }
  
  B_hat = apply( B, mean, MARGIN=c(1,2))
  B_975 = apply( B, quantile, MARGIN=c(1,2),0.975)
  B_025 = apply( B, quantile, MARGIN=c(1,2),0.025)
  colnames(B_hat)=colnames(B_975)=colnames(B_025)=colnames(m$X)
  rownames(B_hat)=rownames(B_975)=rownames(B_025)=colnames(m$Y)
  list(Bsamples = B, Bmean=B_hat, Bq025=B_025,Bq975=B_975)
}



#' Computes and plots the trait-environment relationship of a given CWM trait and a given environmental variable
#'
#' Computes and plots the trait-environment relationship of a given CWM trait and a focal environmental variable. In order to build the response curve, the function builts a dataframe where the focal environmental variable varies along a gradient and the other (non-focal) variables are fixed to their mean (but see FixX parameter for fixing non-focal variables to user-defined values).
#' @param m a model fitted with \code{jtdm_fit}
#' @param indexGradient The name (as specified in the column names of X) of the focal variable.
#' @param indexTrait  The name (as specified in the column names of Y) of the focal trait.
#' @param XFocal Optional. A gradient of the focal variable provided by the user. If provided, the function will used this gradient instead of building a regular one. Default to NULL.
#' @param grid.length The number of points along the gradient of the focal variable. Default to 200.
#' @param FixX Optional. A parameter to specify the value to which non-focal variables are fixed. This can be useful for example if we have some categorical variables (e.g. forest vs meadows) and we want to obtain the partial response curve for a given value of the variable. It has to be a list of the length and names of the columns of X. For example, if the columns of X are "MAT","MAP","Habitat" and we want to fix "Habitat" to 1, then FixX=list(MAT=NULL,MAP=NULL,Habitat=1.). Default to NULL.
#' @param FullPost The type of predictions to be obtain. If FullPost = TRUE, the function returns samples from the predictive distribution. If FullPost="mean", the function computes the posterior distribution of the regression term B\%*\%X). Default to "mean", here FullPost cannot be FALSE.
#' @export
#' @return A list containing:
#'    \item{p}{A plot of the trait-environment relationship.}
#'    \item{predictions}{A data frame containing the predicted trait-environmental relationships including the gradient of the focal environmental variable, mean trait predictions and quantiles (can be useful to code customized plot).}
#'  
#' @examples
#' data(Y)  
#' data(X)  
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10,  
#'         burnin = 100,  
#'         sample = 100)  
#' # SLA-GDD relationship
#' plot = partial_response(m,indexGradient="GDD",indexTrait="SLA")
#' plot$p
#' # SLA-GDD relationship in forest (i.e. when forest=1)
#' plot = partial_response(m,indexGradient="GDD",indexTrait="SLA",
#'                         FixX=list(GDD=NULL,FDD=NULL,forest=1))
#' plot$p
#' @import ggplot2
#' @importFrom stats quantile
#' @importFrom utils globalVariables
partial_response = function(m, indexGradient, indexTrait, XFocal = NULL, grid.length=200,FixX=NULL, FullPost="mean"){
  
  indexGradient = which(colnames(m$X_raw) == indexGradient)
  indexTrait = which(colnames(m$Y) == indexTrait)
  
  if(FullPost==FALSE){stop("FullPost cannon be set to FALSE in partial_response.")}
  if(!is.null(FixX)){if(!identical(names(FixX), colnames(m$X_raw))) {stop("Provide FixX as a list with the same names of X")}}
  if(!is.null(FixX[[names(FixX)[indexGradient]]])) {stop("FixX of the focal environmental variable must be NULL")}
  
  if( !is.null(indexTrait)){if( indexTrait > ncol(m$Y)){stop("index is not available")}}
  if( indexGradient > ncol(m$X)){stop("index is not available")}
  
  data=list(Y=m$Y, X=m$X, K=ncol(m$X), J=ncol(m$Y), n=nrow(m$Y), df= ncol(m$Y), I=diag(ncol(m$Y)),  X_raw = m$X_raw)
  
  #Build XGradient matrix
  if(is.null(XFocal)){
    XGradientFocal=  seq(from=min(data$X_raw[,indexGradient]),to=max(data$X_raw[,indexGradient]),length.out=grid.length)
  }else{XGradientFocal = XFocal
  grid.length = nrow(XFocal)
  }
  
  XGradient_new = matrix(nrow=grid.length,ncol=ncol(data$X_raw))
  
  for(j in 1:ncol(data$X_raw)){
    
    if(j == indexGradient){
      
      XGradient_new[,j] = XGradientFocal
      
    }else{
      
      if(!is.null(FixX[[colnames(data$X_raw)[j]]])){
        XGradient_new[,j] = FixX[[colnames(data$X_raw)[j]]]
      }else{
        XGradient_new[,j] = mean(data$X_raw[,j])
      }
    }
  }
  
  
  colnames(XGradient_new) = colnames(m$X_raw)
  
  PartialPredictions = jtdm_predict(m=m,Xnew=XGradient_new, FullPost=FullPost)
  
  
  table = data.frame(x = XGradientFocal, Predmean = PartialPredictions$PredMean[,indexTrait],
                     Pred975 = PartialPredictions$Predq975[,indexTrait],
                     Pred025 = PartialPredictions$Predq025[,indexTrait])
  
  p = ggplot() + geom_line(data=table, aes(x=x, y=Predmean),col="#00BFC4") +
    geom_ribbon(data=table, aes(x=x,y=Predmean,ymin=Pred025,ymax=Pred975),alpha=0.3) + geom_rug(data=data.frame(x=data$X_raw[,indexGradient]),aes(x=x),sides="b") +
    ggtitle("Partial response curve ") + xlab(colnames(data$X_raw)[indexGradient]) + ylab(colnames(data$Y)[indexTrait]) + theme_minimal()
  
  return(list(p=p,predictions=table))
  
}



#' Get the inferred residual covariance matrix
#'
#' Get the samples from the posterior distribution of the residual covariance matrix, together with the posterior mean and quantiles.
#' @param m a model fitted with \code{jtdm_fit}
#' @export
#' @return A list containing:
#'    \item{Ssamples}{ Sample from the posterior distribution of the residual covariance matrix. It is an array where the first two dimensions are the rows and columns of the matrix, and the third dimensions are the samples from the posterion distribution}
#'   
#'    \item{Smean}{ Posterior mean of the residual covariance matrix.}
#'   
#'    \item{Sq975,Sq025}{ 97.5\% and 0.25\% posterior quantiles of the residual covariance matrix.}
#'
#' @examples
#' data(Y)  
#' data(X) 
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10, 
#'         burnin = 100, 
#'         sample = 100) 
#' # get the inferred residual covariance
#' Sigma =get_sigma(m)
#' @importFrom stats quantile
#' @importFrom coda as.mcmc

get_sigma=function(m){
  
  data=list(Y=m$Y, X=m$X, K=ncol(m$X), J=ncol(m$Y), n=nrow(m$Y), df= ncol(m$Y), I=diag(ncol(m$Y)))
  
  mcmc_param=suppressWarnings(coda::as.mcmc(m$model))
  
  ntot = m$model$sample*length(m$model$mcmc) #samples * n.chains
  
  Sigma = array(dim=c(data$J,data$J,ntot))
  for(i in 1:ntot){
    Sigma[,,i]=matrix(mcmc_param[i,grep("Sigma",colnames(mcmc_param))],ncol=data$J)
  }
  
  Sigma_hat = apply(Sigma, mean, MARGIN=c(1,2))
  Sigma_975 = apply( Sigma, quantile, MARGIN=c(1,2),0.975)
  Sigma_025 = apply( Sigma, quantile, MARGIN=c(1,2),0.025)
  rownames(Sigma_hat)=rownames(Sigma_975)=rownames(Sigma_025)=colnames(Sigma_hat)=colnames(Sigma_975)=colnames(Sigma_025)=colnames(m$Y)
  
  list(Ssamples = Sigma, Smean=Sigma_hat, Sq025=Sigma_025,Sq975=Sigma_975)
  
}



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
#' @details This function is time consuming when \code{FullPost=T}. Consider setting \code{parallel=T} and/or to set \code{mcmc.samples} to a value smaller than the length of the MCMC chains.
#' @examples
#' data(Y)  
#' data(X)  
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 1,
#'         burnin = 10,
#'         sample = 10)  
#' # Compute probability of SLA and LNC to be joint-high at sites in the studies
#' joint = joint_trait_prob(m,indexTrait=c("SLA","LNC"),
#'                          bounds=list(c(mean(Y[,"SLA"]),Inf),c(mean(Y[,"SLA"]),Inf)))
#' @importFrom parallel mclapply detectCores
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats quantile
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



#' Computes partial response curves of joint probabilities
#'
#' Computes the partial responses curves of joint probability of CWM traits as a function of a focal variable. The regions in which joint probabilities are computed are specified by bounds. In order to build the response curve, the function builts a dataframe where the focal variable varies along a gradient and the other (non-focal) variables are fixed to their mean (but see FixX parameter for fixing non-focal variables to user-defined values). Then, uses joint_trait_prob to compute the joint probability in these dataset.
#' @param m A model fitted with \code{jtdm_fit}
#' @param indexTrait  A vector of the names (as specified in the column names of Y) of the two (or more!) traits we want to compute the joint probabilities of.
#' @param indexGradient The name (as specified in the column names of X) of the focal variable.
#' @param bounds The parameter to specify a region in the community-trait space where the function computes the joint probabilities of traits. It is a list of the length of "indexTrait", each element of the list is a vector of length two. The vector represents the inferior and superior bounds of the region for the specified trait. For example, if we consider two traits, bounds=list(c(10,Inf),c(10,Inf)) corresponds to the region in the community-trait space where both traits both take values greater than 10.
#' @param grid.length The number of points along the gradient of the focal variable. Default to 200.
#' @param XFocal Optional. A gradient of the focal variable provided by the user. If provided, the function will used this gradient instead of building a regular one. Default to NULL.
#' @param FixX Optional. A parameter to specify the value to which non-focal variables are fixed. This can be useful for example if we have some categorical variables (e.g. forest vs meadows) and we want to obtain the partial response curve for a given value of the variable. It has to be a list of the length and names of the columns of X. For example, if the columns of X are "MAT","MAP","Habitat" and we want to fix "Habitat" to 1, then FixX=list(MAT=NULL,MAP=NULL,Habitat=1.). Default to NULL.
#' @param FullPost If FullPost = TRUE, the function returns samples from the predictive distribution of joint  probabilities. If FullPost= FALSE, joint probabilities are computed only using the posterior mean of the parameters. FullPost cannot be equal to "mean" here.
#' @param mcmc.samples Optional, default to NULL, only works when FullPost=FALSE. Defines the number of MCMC samples to compute the posterior distribution of joint probabilities. Needs to be between 1 and m$model$sample x length(m$model$mcmc)
#' @param parallel Optional, only works when FullPost = TRUE. When TRUE, the function uses mclapply to parallelise the calculation of the posterior distribution joint probabilities.
#' @details This function is time consuming when \code{FullPost=T}. Consider setting \code{parallel=T} and/or to set \code{mcmc.samples} to a value smaller than the length of the MCMC chains.
#' @export
#' @return A list containing:
#'    \item{GradProbssamples}{Sample from the posterior distribution of the joint probability along the gradient. It is a vector whose lenght is the number of MCMC samples. NULL if FullPost=FALSE. }

#'    \item{GradProbsmean}{Posterior mean of the joint probability along the gradient. }
#'    
#'    \item{GradProbsq975,GradProbsq025}{97.5\% and 0.25\% posterior quantiles of the joint probability along the gradient. NULL if FullPost=FALSE. }
#'    
#'    \item{gradient}{The gradient of the focal variable built by the function.}
#' 
#' @examples
#' data(Y)  
#' data(X)  
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 1,  
#'         burnin = 10,  
#'         sample = 10)  
#' # Compute probability of SLA and LNC to be joint-high at sites in the studies
#'
#' # Compute the joint probability of SLA and LNC 
#'   #to be joint-high along the GDD gradient
#' joint = joint_trait_prob_gradient(m,indexTrait=c("SLA","LNC"), 
#'                                   indexGradient="GDD",
#'                                   bounds=list(c(mean(Y[,"SLA"]),Inf),c(mean(Y[,"SLA"]),Inf)))
#' # Compute the joint probability of SLA and LNC to be joint-high along the
#' # GDD gradient when forest = 1 (i.e. in forests) 
#' joint = joint_trait_prob_gradient(m,indexTrait=c("SLA","LNC"),
#'                                   indexGradient="GDD",
#'                                   bounds=list(c(mean(Y[,"SLA"]),Inf),c(mean(Y[,"SLA"]),Inf)),
#'                                   FixX=list(GDD=NULL,FDD=NULL,forest=1))
#'                      
#' @importFrom stats quantile 
#' @importFrom parallel mclapply detectCores

joint_trait_prob_gradient = function(m, indexTrait, indexGradient,bounds, grid.length=200, XFocal=NULL, FixX=NULL, FullPost=T,mcmc.samples=NULL,parallel=FALSE){
  
  indexTrait = sapply(indexTrait,function(x){which(colnames(m$Y) %in% x )})
  indexGradient = which(colnames(m$X_raw) == indexGradient)
  
  ntot = m$model$sample*length(m$model$mcmc) #samples * n.chains
  if(is.null(mcmc.samples)){mcmc.samples=ntot}
  if(mcmc.samples>ntot){stop("You need to provide a number of mcmc samples lower than the length of the chain given by m$model$sample*length(m$model$mcmc)")}
  if(FullPost=="mean"){stop("Fullstop cannot be set to mean here.")}
  if(length(indexTrait) != length(bounds)){stop("index and bounds have different lengths!!")}
  if(length(indexGradient)>1) {stop("indexGradient has to be a numerical value")}
  if(!is.null(FixX)){if(!identical(names(FixX), colnames(m$X_raw))) {stop("Provide FixX as a list with the same names of X")}}
  if(!is.null(FixX[[names(FixX)[indexGradient]]])) {stop("FixX of the focal environmental variable must be NULL")}
  
  data=list(Y=m$Y, X=m$X, K=ncol(m$X), J=ncol(m$Y), n=nrow(m$Y), df= ncol(m$Y), I=diag(ncol(m$Y)), X_raw = m$X_raw)
  
  
  #Build XGradient matrix
  if(is.null(XFocal)){
    XGradientFocal=  seq(from=min(data$X_raw[,indexGradient]),to=max(data$X_raw[,indexGradient]),length.out=grid.length)
  }else{XGradientFocal = XFocal
  grid.length = nrow(XFocal)
  }
  
  XGradient = matrix(nrow=grid.length,ncol=ncol(data$X_raw))
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
                                 Xnew=XGradient, bounds=bounds, FullPost=T, 
                                 mcmc.samples = mcmc.samples, parallel = parallel)$PROBsamples
    GradProbs_hat =  apply(GradProbs,mean,MARGIN=1)
    GradProbs_975 = apply( GradProbs, quantile, MARGIN=1,0.975)
    GradProbs_025 = apply( GradProbs, quantile, MARGIN=1,0.025)
    
  }else{
    GradProbs_hat = joint_trait_prob(m, indexTrait=colnames(m$Y)[indexTrait],
                                     Xnew=XGradient, bounds=bounds, FullPost=F)$PROBmean
    GradProbs_975 = GradProbs_025 = NULL
  }
  
  
  GradProbs_hat =  apply(GradProbs,mean,MARGIN=1)
  GradProbs_975 = apply( GradProbs, quantile, MARGIN=1,0.975)
  GradProbs_025 = apply( GradProbs, quantile, MARGIN=1,0.025)
  list(GradProbssamples = GradProbs, GradProbsmean=GradProbs_hat,
       GradProbsq025=GradProbs_025,GradProbsq975=GradProbs_975,gradient=XGradientFocal)
  
}



#' jtdm.
#'
#' Package to fit a Join Trait Distribution Model and to analyse its result to understand and predict the community-level strategy. See Poggiato et al. In prep.
#'
#'
#' @docType package
#'
#'
#' @name jtdm
NULL

#' Site x environmental covariates dateset
#'
#' Includes the Growing Degree Days (GDD) during the growing season and Freezing Degree Days (FDD) during the growing season averaged over the period 1989-2019
#'
#' @docType data
#'
#' @usage data(X)
#'
#' @format A matrix
#'
#' @usage data(X)
#'
#' @name X
#'
#' @author Orchamp consortium
#'
#' @keywords datasets
#'
#' @examples
#' data(X)
NULL

#' Site x CWM traits dataset
#'
#' A site x CWM traits dataset computed using pinpoint abundances of plants and species mean
#  traits
#'
#' @docType data
#'
#' @usage data(Y)
#'
#' @format A matrix
#'
#' @author Orchamp Consortium
#'
#' @keywords datasets
#'
#' @name Y
#'
#' @examples
#' data(Y)
NULL



#' Partial response curve of the pairwise most suitable community-level strategy and of the pairwise envelop of possible community-level strategy
#'
#' Partial response curve of the pairwise most suitable community-level strategy and of the pairwise envelop of possible community-level strategy. In order to build the response curve, the function builts a dataframe where the focal variable varies along a gradient and the other (non-focal) variables are fixed to their mean (but see FixX parameter for fixing non-focal variables to user-defined values). The chosen traits are specified in indexTrait. Then uses the jtdm_predict function to compute the most suitable community-level strategy and the residual covariance matrix to build the envelop of possible CLS.
#' @param m a model fitted with \code{jtdm_fit}
#' @param indexTrait  A vector of the two names (as specified in the column names of Y) containing the two (or more!) traits we want to compute the community level strategy of.
#' @param indexGradient The name (as specified in the column names of X) of the focal variable.
#' @param FullPost   If FullPost = TRUE, the function returns samples from the predictive distribution of joint probabilities. If FullPost= FALSE, joint probabilities are computed only using the posterior mean of the parameters.
#' @param grid.length The number of points along the gradient of the focal variable. Default to 20 (which ensures a fair visualisation).
#' @param FixX Optional. A parameter to specify the value to which non-focal variables are fixed. This can be useful for example if we have some categorical variables (e.g. forest vs meadows) and we want to obtain the partial response curve for a given value of the variable. It has to be a list of the length and names of the columns of X. For example, if the columns of X are "MAT","MAP","Habitat" and we want to fix "Habitat" to 1, then FixX=list(MAT=NULL,MAP=NULL,Habitat=1.). Default to NULL.
#' @param confL The confidence level of the confidence ellipse (i.e. of the envelop of possible community-level strategies). Default is 0.95.
#' @return Plot of the partial response curve of the pairwise most suitable community-level strategy and of the pairwise envelop of possible community-level strategy
#' @export
#' @examples
#' data(Y)  
#' data(X)  
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10,  
#'         burnin = 100,  
#'         sample = 100)  
#'
#' # plot the pairwise SLA-LNC partial response curve along the GDD gradient
#' ellipse_plot(m,indexTrait = c("SLA","LNC"),indexGradient="GDD")
#' #  plot the pairwise SLA-LNC partial response curve along the GDD gradient
#' #  in forest (i.e. when forest=1)
#' ellipse_plot(m,indexTrait = c("SLA","LNC"),indexGradient="GDD",
#'              FixX=list(GDD=NULL,FDD=NULL,forest=1))
#' @importFrom stats qchisq
#' @import ggplot2
#' @importFrom ggforce geom_ellipse
ellipse_plot = function(m,indexGradient,indexTrait,FullPost=F, grid.length=20, FixX=NULL, confL= 0.95){
  
  indexGradient = which(colnames(m$X_raw) == indexGradient)
  indexTrait = sapply(indexTrait,function(x){which(colnames(m$Y) %in% x )})
  
  if(!is.null(FixX)){if(!identical(names(FixX), colnames(m$X_raw))) {stop("Provide FixX as a list with the same names of X")}}
  if(!is.null(FixX[[names(FixX)[indexGradient]]])) {stop("FixX of the focal environmental variable must be NULL")}
  
  
  if(length(indexGradient)>1) {stop("indexGradient has to be a scalar")}
  if(length(indexTrait)!=2) {stop("indexTrait has to be a vector of length 2")}
  
  data=list(Y=m$Y, X=m$X, K=ncol(m$X), J=ncol(m$Y), n=nrow(m$Y), df= ncol(m$Y), I=diag(ncol(m$Y)),  X_raw = m$X_raw)
  
  Sigma=get_sigma(m)
  
  #create gradient of one variable
  XGradientFocal =  seq(from=min(data$X_raw[,indexGradient]),to=max(data$X_raw[,indexGradient]),length.out=grid.length)
  
  
  XGradient_new = matrix(nrow=grid.length,ncol=ncol(data$X_raw))
  
  for(j in 1:ncol(data$X_raw)){
    
    if(j == indexGradient){
      
      XGradient_new[,j] = XGradientFocal
      
    }else{
      
      if(!is.null(FixX[[colnames(data$X_raw)[j]]])){
        XGradient_new[,j] = FixX[[colnames(data$X_raw)[j]]]
      }else{
        XGradient_new[,j] = mean(data$X_raw[,j])
      }
    }
    
    
  }
  
  colnames(XGradient_new) = colnames(data$X_raw)
  
  
  
  #plot for XGradient_new
  if(FullPost){
    #If we want to extract the posterior distribution
    prediction = jtdm_predict(m=m, Xnew = XGradient_new, validation = F, FullPost=T)
    table=data.frame(X = XGradientFocal, prediction$PredMean[,indexTrait])
  }else{
    prediction =jtdm_predict(m=m, Xnew=XGradient_new, validation = F,FullPost=F)
    table=data.frame(X = XGradientFocal, prediction$PredMean[,indexTrait])
  }
  
  eigen1 = eigen(Sigma$Smean[indexTrait,indexTrait])$value[1]
  eigen2 = eigen(Sigma$Smean[indexTrait,indexTrait])$value[2]
  angle = atan(eigen(Sigma$Smean[indexTrait,indexTrait])$vectors[2,1]/eigen(Sigma$Smean[indexTrait,indexTrait])$vectors[1,1])
  # How to define ellipses parameters?
  # http://www.cs.utah.edu/~tch/CS6640F2020/resources/How%20to%20draw%20a%20covariance%20error%20ellipse.pdf
  # or http://jmlilly.net/course/v3/pdfs/thevarianceellipse.pdf
  # https://www.xarg.org/2018/04/how-to-plot-a-covariance-error-ellipse/
  
  p=ggplot(data = table,
           mapping = aes_string(x = colnames(data$Y)[indexTrait[1]],
                                y = colnames(data$Y)[indexTrait[2]])) +
    geom_ellipse(
      aes_string(
        x0 = colnames(data$Y)[indexTrait[1]],
        y0 = colnames(data$Y)[indexTrait[2]],
        fill = "X",
        a = sqrt(eigen1 * qchisq(confL,df=2)),
        b = sqrt(eigen2 * qchisq(confL,df=2)),
        angle = angle
      ),
      col = "NA",
      alpha = 0.1
    ) +
    geom_point(shape = 21, size = 3, aes_string(fill = "X")) +
    labs(
      x = colnames(data$Y)[indexTrait[1]],
      y = colnames(data$Y)[indexTrait[2]],
      color = colnames(data$X_raw)[indexGradient]
    ) +
    geom_ellipse(
      aes_string(
        x0 = colnames(data$Y)[indexTrait[1]],
        y0 = colnames(data$Y)[indexTrait[2]],
        a = sqrt(eigen1 * qchisq(confL,df=2)),
        b = sqrt(eigen2 * qchisq(confL,df=2)),
        angle = angle
      ),
      col = "black",
      size = 0.25
    ) +
    guides(fill = guide_legend(title = colnames(data$X_raw)[indexGradient])) +
    scale_fill_viridis_c() +
    theme_classic()
  
  return(p)
  
  
}




#' K-fold cross validation predictions and goodness of fit metrics
#'
#' Run K-fold cross validation predictions of the model m on a specified dataset.
#' @param m a model fitted with \code{jtdm_fit}
#' @param K  The number of folds of the K-fold cross validation
#' @param partition A partition of the dataset specified by the user. It is a vector (whose length are the number of sites), where each element specifies the fold index of the site.
#' @param adapt,burnin,sample,n.chains  Parameters of the MCMC sampler. See \code{?run.jags} for details
#' @export
#' @return A list containing:
#'    \item{Pred}{Sample from the posterior predictive distribution in cross validation. It is an array where the first dimension is the number of sites in Xnew, the second is the number of traits modelled and the third the number of MCMC samples. NULL if FullPost=FALSE. }
#'    
#'    \item{PredMean}{Posterior mean of posterior predictive distribution in cross validation. }
#'    
#'    \item{Predq975,Predq025}{97.5\% and 0.25\% posterior quantiles of the posterior predictive distribution in cross validation. NULL if FullPost=FALSE. }
#'    
#'    \item{R2}{R squared of predictions in cross validation. }
#'    
#'    \item{RMSE}{Root square mean error between  squared of predictions in cross validation.}
#'
#' @examples
#' data(Y)  
#' data(X)  
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 1,  
#'         burnin = 10,  
#'         sample = 10)  
#' # Run 3-fold cross validation on m
#' pred = jtdmCV(m, K=5, adapt = 1, burnin = 10, sample = 10)
#' @importFrom stats quantile

jtdmCV = function(m, K = 5,
                  adapt = 200,
                  burnin = 500,
                  sample = 500,
                  n.chains=2, partition=NULL){
  
  data=list(Y=m$Y, X=m$X, K=ncol(m$X), J=ncol(m$Y), n=nrow(m$Y), df= ncol(m$Y), I=diag(ncol(m$Y)),  X_raw = m$X_raw)
  
  # Create partition (if needed)
  if(!is.null(partition)){if(length(partition) != data$n) stop("partition must be a vector of length n (the number of sites")
  }else{index <- sample(1:K,size=data$n,replace=TRUE,prob=rep(data$n/K,K))}
  
  
  ntot = sample*n.chains
  preds = array(dim=c(data$n,data$J,ntot))
  RMSE = R2 = matrix(nrow=K,ncol=data$J)
  
  # CV loop
  for(i in 1:K){
    train = which(index != i)
    test = which(index == i)
    Y=data$Y[train,]
    X_raw=data$X_raw[train,]
    
    m = jtdm_fit(Y=Y, X=X_raw,
                 formula=m$formula,
                 adapt = adapt,
                 burnin = burnin,
                 sample = sample,
                 n.chains = n.chains,
                 monitor = c('B','Sigma'))
    
    prediction = jtdm_predict(m=m,Xnew=data$X_raw[test,],Ynew=data$Y[test,],validation = T,FullPost = T)
    preds[test,,] = prediction$Pred
    
    RMSE[i,] = prediction$RMSE
    R2[i,] = prediction$R2
    
    cat("Fold ", i, " out of ", K,"\n")
    
  }
  
  meanPred = apply(preds,mean,MARGIN = c(1,2))
  Pred975 = apply(preds, quantile, MARGIN=c(1,2),0.975)
  Pred025 = apply(preds, quantile, MARGIN=c(1,2),0.025)
  
  R2.cv = apply(R2,mean,MARGIN = 2)
  names(R2.cv) = colnames(data$Y)
  
  RMSE.cv = apply(RMSE,mean,MARGIN = 2)
  names(RMSE.cv) = colnames(data$Y)
  
  list(Pred = preds, PredMean = meanPred, Predq025 = Pred025,
       Predq975 = Pred975, R2 = R2.cv,
       RMSE = apply(RMSE,mean,MARGIN = 2), partition = index)
  
}
