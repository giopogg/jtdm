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

  list(Pred = preds, PredMean = meanPred, Predq025 = Pred025, Predq975 = Pred975, R2 = R2.cv, RMSE = apply(RMSE,mean,MARGIN = 2), partition = index)

}
