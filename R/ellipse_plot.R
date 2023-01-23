#' Partial response curve of the pairwise most suitable community-level strategy and of the pairwise envelop of possible community-level strategy
#'
#' Partial response curve of the pairwise most suitable community-level strategy and of the pairwise envelop of possible community-level strategy. In order to build the response curve, the function builds a dataframe where the focal variable varies along a gradient and the other (non-focal) variables are fixed to their mean (but see FixX parameter for fixing non-focal variables to user-defined values). The chosen traits are specified in indexTrait. Then uses the jtdm_predict function to compute the most suitable community-level strategy and the residual covariance matrix to build the envelop of possible CWM combinations.
#' @param m a model fitted with \code{jtdm_fit}
#' @param indexTrait  A vector of the two names (as specified in the column names of Y) containing the two (or more!) traits we want to compute the community level strategy of.
#' @param indexGradient The name (as specified in the column names of X) of the focal variable.
#' @param FullPost If FullPost = TRUE, the function returns samples from the predictive distribution of joint probabilities. If FullPost= FALSE, joint probabilities are computed only using the posterior mean of the parameters.
#' @param grid.length The number of points along the gradient of the focal variable. Default to 20 (which ensures a fair visualization).
#' @param FixX Optional. A parameter to specify the value to which non-focal variables are fixed. This can be useful for example if we have some categorical variables (e.g. forest vs meadows) and we want to obtain the partial response curve for a given value of the variable. It has to be a list of the length and names of the columns of X. For example, if the columns of X are "MAT","MAP","Habitat" and we want to fix "Habitat" to 1, then FixX=list(MAT=NULL,MAP=NULL,Habitat=1.). Default to NULL.
#' @param confL The confidence level of the confidence ellipse (i.e. of the envelop of possible community-level strategies). Default is 0.95.
#' @return Plot of the partial response curve of the pairwise most suitable community-level strategy and of the pairwise envelop of possible community-level strategy
#' @export
#' @examples
#' data(Y)  
#' data(X)  
#' # Short MCMC to obtain a fast example: results are unreliable !
#' m = jtdm_fit(Y=Y, X=X, formula=as.formula("~GDD+FDD+forest"),  sample = 1000)  
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
ellipse_plot = function(m, indexGradient, indexTrait, FullPost = FALSE, grid.length = 20, FixX = NULL, confL = 0.95){

  if(!inherits(m, "jtdm_fit")) stop("m is not an object of class jtdm_fit")
  
  indexGradient = which(colnames(m$X_raw) == indexGradient)
  indexTrait = sapply(indexTrait,function(x){which(colnames(m$Y) %in% x )})

  if(!is.null(FixX)){if(!identical(names(FixX), colnames(m$X_raw))) {
    stop("Provide FixX as a list with the same names of X")}}
  if(!is.null(FixX[[names(FixX)[indexGradient]]])) {
    stop("FixX of the focal environmental variable must be NULL")}


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
    prediction = jtdm_predict(m=m, Xnew = XGradient_new, validation = FALSE, FullPost = TRUE)
    table=data.frame(X = XGradientFocal, prediction$PredMean[,indexTrait])
  }else{
    prediction =jtdm_predict(m=m, Xnew=XGradient_new, validation = FALSE,FullPost = FALSE)
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
