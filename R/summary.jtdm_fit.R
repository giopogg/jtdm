#' Prints the summary of a fitted jtdm
#' 
#' Prints the summary of a fitted jtdm
#' @param object a model fitted with \code{jtdm_fit}
#' @param ... 	additional arguments
#' @return A printed summary of the fitted jtdm
#' @author Giovanni Poggiato
#' @examples
#' data(Y)  
#' data(X)  
#' m = jtdm_fit(Y=Y, X=X, 
#'              formula=as.formula("~GDD+FDD+forest"), sample = 1000)  
#' summary(m)
#' @method summary jtdm_fit
#' @export
#' 

summary.jtdm_fit = function(object, ...){
  
  m = object
  if(!inherits(m, "jtdm_fit")) stop("object is not an object of class SDMfit" )
  cat("================================================================== \n")
  
  model = paste0("JTDM to model ", paste0(colnames(m$Y),collapse = ", "), "\n",
                 "as a function of ", sub("~","", m$formula)[2],".\n",
                 "Likelihood of the model: ", m$model$log.lik, " \n")
  cat(model)
  cat("================================================================== \n")
  cat("* Useful S3 methods and functions to get and show inferred parameters: \n")
  cat("plot(), getB(), get_sigma() \n")
  cat("* Useful functions to predict marginally and jointly: \n")
  cat("jtdm_predict(), jtdmCV(), joint_trait_prob() \n")
  cat("* Function to plot partial response curves of the envelope of possible: \n")
  cat("ellipse_plot()\n")
  cat("* Function to compute partial response curves of joint probabilities: \n")
  cat("joint_trait_prob_gradient() \n")
  cat("================================================================== \n")

  #Just to fix pkgdown
  invisible(object)
  
}