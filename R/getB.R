#' Get the inferred regression coefficients
#'
#' Get the samples from the posterior distribution of the regression coefficient matrix B, together with the posterior mean and quantiles. The regression coefficient matrix B is a matrix where the number of rows is defined by the number of traits that are modeled, and the number of columns is the number of columns of the matrix m$X (the number of explanatory variables after transformation via formula)
#' @param m a model fitted with \code{jtdm_fit}
#' @export
#' @return A list containing:
#'    \item{Bsamples}{Sample from the posterior distribution of the regression coefficient matrix. It is an array where the first dimension is the number of traits, the second the number of columns in m$X (the number of variables after transformation via formula) and the third the number of MCMC samples.}
#'    \item{Bmean}{Posterior mean of the regression coefficient matrix.}
#'    \item{Bq975,Bq025}{97.5\% and 0.25\% posterior quantiles of the regression coefficient matrix.}
#' @examples
#' data(Y)  
#' data(X)
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"), sample = 1000) 
#' # get the inferred regression coefficients
#' B=getB(m)
#' @importFrom stats quantile

getB=function(m){

  if(!inherits(m, "jtdm_fit")) stop("m is not an object of class jtdm_fit")
  data=list(Y = m$Y, X = m$X, K = ncol(m$X), J = ncol(m$Y), n = nrow(m$Y), df= ncol(m$Y), I=diag(ncol(m$Y)))

  B = m$model$B

  B_hat = apply( B, mean, MARGIN = c(1,2))
  B_975 = apply( B, quantile, MARGIN = c(1,2), 0.975)
  B_025 = apply( B, quantile, MARGIN = c(1,2), 0.025)
  colnames(B_hat) = colnames(B_975) = colnames(B_025) = colnames(m$X)
  rownames(B_hat) = rownames(B_975) = rownames(B_025) = colnames(m$Y)
  list(Bsamples = B, Bmean=B_hat, Bq025=B_025,Bq975=B_975)
}
