#' Computes and plots the trait-environment relationship of a given CWM trait and a given environmental variable
#'
#' Computes and plots the trait-environment relationship of a given CWM trait and a focal environmental variable. In order to build the response curve, the function builds a dataframe where the focal environmental variable varies along a gradient and the other (non-focal) variables are fixed to their mean (but see FixX parameter for fixing non-focal variables to user-defined values).
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
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"), sample = 1000)  
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
partial_response = function(m, indexGradient, indexTrait, XFocal = NULL,
                            grid.length = 200,FixX = NULL, FullPost = "mean"){
  
  if(!inherits(m, "jtdm_fit")) stop("m is not an object of class jtdm_fit")
  
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
    XGradientFocal = seq(from = min(data$X_raw[,indexGradient]),
                         to = max(data$X_raw[,indexGradient]),
                         length.out = grid.length)
  }else{
    XGradientFocal = XFocal
    grid.length = nrow(XFocal)
  }
  
  XGradient_new = matrix(nrow = grid.length, ncol = ncol(data$X_raw))
  
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
