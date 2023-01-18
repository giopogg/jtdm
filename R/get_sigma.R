#' Get the inferred residual covariance matrix
#'
#' Get the samples from the posterior distribution of the residual covariance matrix, together with the posterior mean and quantiles.
#' @param m a model fitted with \code{jtdm_fit}
#' @export
#' @return A list containing:
#'    \item{Ssamples}{ Sample from the posterior distribution of the residual covariance matrix. It is an array where the first two dimensions are the rows and columns of the matrix, and the third dimensions are the samples from the posterior distribution}
#'   
#'    \item{Smean}{ Posterior mean of the residual covariance matrix.}
#'   
#'    \item{Sq975,Sq025}{ 97.5\% and 0.25\% posterior quantiles of the residual covariance matrix.}
#'
#' @examples
#' data(Y)  
#' data(X) 
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"), sample = 1000) 
#' # get the inferred residual covariance
#' Sigma =get_sigma(m)
#' @importFrom stats quantile

get_sigma=function(m){
  
  if(!inherits(m, "jtdm_fit")) stop("m is not an object of class jtdm_fit")
  
  data=list(Y = m$Y, X = m$X, K = ncol(m$X), J = ncol(m$Y), n = nrow(m$Y), df = ncol(m$Y), I = diag(ncol(m$Y)))

  Sigma = m$model$Sigma
  
  Sigma_hat = apply(Sigma, mean, MARGIN = c(1,2))
  Sigma_975 = apply( Sigma, quantile, MARGIN = c(1,2),0.975)
  Sigma_025 = apply( Sigma, quantile, MARGIN = c(1,2),0.025)
  rownames(Sigma_hat) = rownames(Sigma_975) = rownames(Sigma_025) = colnames(Sigma_hat) = colnames(Sigma_975) = colnames(Sigma_025) = colnames(m$Y)

  list(Ssamples = Sigma, Smean=Sigma_hat, Sq025=Sigma_025,Sq975=Sigma_975)

}
