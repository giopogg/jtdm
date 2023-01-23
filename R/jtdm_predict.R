#' Predict method for joint trait distribution model
#'
#' Obtains predictions from a fitted joint trait distribution model and optionally computes their R squared and root mean square error (RMSE)
#' @param m a model fitted with \code{jtdm_fit}
#' @param Xnew optionally, a data frame in which to look for variables with which to predict. If omitted, the fitted linear predictors are used
#' @param validation  boolean parameter to decide whether we want to compute goodness of fit measures. If true, then Ynew is needed.
#' @param Ynew   Optional. The observed response variables at sites specified in Xnew. It is used to compute goodness of fit metrics when validation= T.
#' @param FullPost   The type of predictions to be obtain. If FullPost = TRUE, the function returns samples from the predictive distribution, the credible intervals are thus the predictive credible interval. If FullPost="mean", the function computes the posterior distribution of the regression term \eqn{BXnew}), i.e., classical credible intervals. If FullPost=FALSE, the function only returns the posterior mean of the regression term (\eqn{BmeanXnew}), i.e., no credible intervals. 
#' @details To obtain a full assessment of the posterior distribution, the function should be ran with FullPost=TRUE, although this can be time consuming. FullPost="mean" is used to compute partial response curves, while FullPost=FALSE is used to compute goodness of fit metrics.
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
#' m = jtdm_fit(Y = Y, X = X, formula=as.formula("~GDD+FDD+forest"), sample = 1000)
#' # marginal predictions of traits in the sites of X
#' pred = jtdm_predict(m)
#' @importFrom stats model.frame model.matrix quantile cor 
#' @importFrom mniw rmNorm
jtdm_predict = function(m = m, Xnew = NULL, Ynew = NULL, validation = FALSE, FullPost = "mean"){

  if(!inherits(m, "jtdm_fit")) stop("m is not an object of class jtdm_fit")
  
  data=list(Y=m$Y, X=m$X, K=ncol(m$X), J=ncol(m$Y), n=nrow(m$Y), df= ncol(m$Y), I=diag(ncol(m$Y)),  X_raw = m$X_raw)

  if(is.null(Xnew)) Xnew=m$X_raw
  if(is.null(dim(Xnew))){Xnew=t(as.matrix(Xnew))}
  if(validation == TRUE & is.null(Ynew)) stop(" if validation = T, you need to provide Ynew!")
  if(ncol(Xnew) != ncol(data$X_raw)) stop("The number of columns of X and Xnew differ")
  if(!is.null(Ynew)){
    if(nrow(Xnew) != nrow(Ynew)) stop("The number of line of Xnew and Ynew differ")
    if(ncol(Ynew) != ncol(data$Y)) stop("The number of columns of Y and Ynew differ")
  }
  if(is.null(FullPost)) stop("please tell whether you want predictions or fitted")
  if(!identical(colnames(Xnew),colnames(m$X_raw))) stop("Provide same column names and same order of the colums in Xnew!")

  ###### trasform Xnew with formula
  Xnew_raw = Xnew
  Xnew = model.frame(m$mt,as.data.frame(Xnew))
  Xnew = model.matrix(m$mt,Xnew)

  ### Compute predictions
  if(FullPost != FALSE){

    B = m$model$B
    Sigma = m$model$Sigma
    ntot = m$model$sample*length(m$model$mcmc) #samples * n.chains

    Predictions = array(dim = c(nrow(Xnew), ncol(data$Y), dim(B)[3]))
    
    for(i in 1: dim(B)[3]){
      if(FullPost == TRUE){
        temp_pred = rmNorm(nrow(Xnew), Xnew %*% t(B[,,1]), Sigma[,,i])
        rownames(temp_pred) = rownames(Xnew)
        Predictions[,,i] = temp_pred
      }else{
        Predictions[,,i] = Xnew %*% t(B[,,i])
      }
    }
    
    meanPred = apply(Predictions, mean, MARGIN = c(1,2))
    Pred975 = apply(Predictions, quantile, MARGIN = c(1,2),0.975)
    Pred025 = apply(Predictions, quantile, MARGIN = c(1,2),0.025)
    colnames(meanPred) = colnames(Pred975) = colnames(Pred025) = colnames(data$Y)
    if(!is.null(rownames(Xnew))) rownames(meanPred) = rownames(Pred975) = rownames(Pred025) = rownames(Xnew)

  }else{ #If we only want to compoute the mean
    
    meanPred = as.matrix(Xnew) %*% t(as.matrix(getB(m)$Bmean))

    colnames(meanPred) = colnames(data$Y)
    if(!is.null(rownames(Xnew))) rownames(meanPred) = rownames(Xnew)

    Pred025 = Pred975 = NULL
    Predictions = NULL
  }

  R2_mod = RMSE_mod = NULL
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
