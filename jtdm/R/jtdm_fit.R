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
