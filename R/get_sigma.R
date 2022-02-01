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
#' }
#' @examples
#' data(Y)  
#' data(X) 
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  adapt = 10, 
#'         burnin = 100, 
#'         sample = 100) 
#' # get the inferred residual covariance
#' Sigma =get_sigma(m)


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
