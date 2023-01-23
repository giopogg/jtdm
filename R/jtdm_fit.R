#' Fitting joint trait distribution models
#'
#' jtdm_fit is used to fit a Joint trait distribution model. Requires the response variable Y (the sites x traits matrix) and the explanatory variables X.This function samples from the posterior distribution of the parameters, which has been analytically determined. Therefore, there is no need for classical MCMC convergence checks. 
#' @param Y The sites x traits matrix containing community (weighted) means of each trait at each site.
#' @param X The design matrix, i.e. sites x predictor matrix containing the value of each explanatory variable (e.g. the environmental conditions) at each site.
#' @param formula An object of class "formula" (or one that can be coerced to that class): a symbolic description of the model to be fitted. The details of model specification are given under 'Details'.
#' @param sample Number of samples from the posterior distribution. Since we sample from the exact posterior distribution, the number of samples is relative lower than MCMC samplers. As a rule of thumb, 1000 samples should provide correct inference.
#' @export
#' @details A formula has an implied intercept term. To remove this use either y ~ x - 1 or y ~ 0 + x. See formula for more details of allowed formulae.
#' @return A list containing:
#'    \item{model}{ An object of class "jtdm_fit", containing the samples from the posterior distribution of the regression coefficients (B) and residual covariance matrix (Sigma), together with the likelihood of the model.}
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
#' m = jtdm_fit(Y = Y, X = X, formula = as.formula("~GDD+FDD+forest"), sample = 1000)  
#' @importFrom stats model.frame model.matrix rWishart coef
#' @importFrom mniw rMT riwish

jtdm_fit = function(Y, X,
                formula,
                sample = 1000
                ){

  if(nrow(Y) != nrow(X)) stop("The number of lines of X and Y do not coincide")
  
  X_raw = X
  X=model.frame(formula,as.data.frame(X))
  mt <- attr(X, "terms")
  X <- model.matrix(mt, X)

  # data preparation
  data=list(Y=Y, X=X, K=ncol(X), J=ncol(Y), n=nrow(Y), df= ncol(Y), I=diag(ncol(Y)))

  # Define prior hyperparameters
  n = data$n
  q = data$K
  p = nu = data$df
  B_0 = matrix(0, nrow = p, ncol = q)
  D = diag(q)*10^4
  Q = diag(ncol(Y)) # Probably needs to play
   
  #########################################################################################################
  ### Sample from the conjugate posterior (see Rowe 2002)
  
  # Posterior hyperparameters of B
  df_post = n + nu - p - 1
  B_bar = ( t(Y) %*% X + B_0%*%solve(D) ) %*% solve( solve(D) + t(X) %*% X)
  G = Q + t(Y)%*%Y + B_0 %*% solve(D) %*% t(B_0) - (t(Y) %*% X + B_0 %*% solve(D)) %*% solve(solve(D) + t(X) %*% X) %*% t(t(Y) %*% X + B_0 %*% solve(D))
  
  # Sample B
  B = rMT(n = sample,
      Lambda = B_bar,
      SigmaC = solve(df_post * (solve(D) + t(X) %*% X)),
      SigmaR = G,
      nu = df_post)
  
  # Posterior hyperparameters of Sigma
  nu_post = n + nu
  Sigma = riwish(sample,
                nu = n + nu,
                Psi = G)
     
  #########################################################################################################
  ### Compute the likelihood (using posterior means as estimates of the parameters)
  
  Sigma_bar =  G/(n + nu - 2*p -2)

  log.lik = -(n*p/2)*log(2*pi) - (n/2)*log(det(Sigma_bar)) -1/2*sum(diag((Y-X%*%t(B_bar))%*%solve(Sigma_bar)%*%t((Y-X%*%t(B_bar)))))

  fitted_jtdm = list(model= list(B = B, Sigma = Sigma, log.lik = log.lik),
                     Y=Y, X=X, X_raw = X_raw, formula=formula,mt=mt)

  class(fitted_jtdm) = "jtdm_fit"
  
  fitted_jtdm
  
}
