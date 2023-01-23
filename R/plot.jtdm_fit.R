#' Plots the parameters of a fitted jtdm
#' 
#' Plots the regression coefficients and covariance matrix of a fitted jtdm
#' @param x a model fitted with \code{jtdm_fit}
#' @param ... 	additional arguments
#' @return A plot of the regression coefficients and covariance matrix of the fitted model
#' @author Giovanni Poggiato
#' @examples
#' data(Y)  
#' data(X)  
#' m = jtdm_fit(Y=Y, X=X, 
#'              formula=as.formula("~GDD+FDD+forest"), sample = 1000)  
#' plot(m)
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @importFrom reshape2 melt
#' @importFrom  grDevices devAskNewPage
#' @importFrom stats cov2cor
#' @method plot jtdm_fit
#' @export
#' 

plot.jtdm_fit = function(x, ...){
  
  m = x
  
  if(!inherits(m, "jtdm_fit")) stop("x is not an object of class SDMfit" )
  
  B = getB(m)
  
  table = data.frame(Bm = as.vector(B$Bmean),
                     B97=as.vector(B$Bq975), 
                     B02 = as.vector(B$Bq025),
                     trait = rep(rownames(B$Bmean),ncol(m$X)),
                     predictor = rep(colnames(B$Bmean), each = nrow(B$Bmean))
                     
  )
  
  table[,"significant"] = ifelse(sign(table$B97)==sign((table$B02)),"yes","no")
  

  p_B = ggplot(data = table, aes(x = Bm, y = predictor, color = significant)) +
    geom_point(size=2) +
    geom_errorbarh(aes(xmax = B97, xmin = B02, height = 0)) +
    geom_vline(xintercept=0,linetype="dashed") +
    facet_grid(trait~.) +
    theme_minimal()+ 
    theme(axis.title.y = element_blank(),axis.title.x=element_blank())
  
  Sigma = get_sigma(m = m)
  
  Sigma_sign = ifelse(sign(Sigma$Sq025)==sign(Sigma$Sq975),1,0)
  Sigma_plot = cov2cor(Sigma$Smean) * Sigma_sign
  colnames(Sigma_plot) = rownames(Sigma_plot) = colnames(Y)
  
  # Reproduce code from ggcorrplot
  
  corr = base::round(x = Sigma_plot, digits = 2)
  corr[which(lower.tri(corr, diag = FALSE) == FALSE)] = NA
  corr <- reshape2::melt(corr, na.rm = TRUE)
  colnames(corr) <- c("Var1", "Var2", "value")
  label <- corr[, "value"]
  
  p_S <- ggplot2::ggplot(data = corr, mapping = ggplot2::aes_string(x = "Var1", 
                                                                    y = "Var2", fill = "value")) +
    geom_tile(color = "gray") +
    scale_fill_gradient2(low = "blue", high = "red", 
                         mid = "white", midpoint = 0, limit = c(-1, 1), space = "Lab", 
                         name = "Correlation") +
    geom_text(mapping = aes_string(x = "Var1", y = "Var2"), label = label, color = "black", size = 5) + theme_minimal() +
    theme(axis.title.y = element_blank(),axis.title.x=element_blank())
  
  
  if(ncol(m$Y) > 6){
    devAskNewPage(TRUE)
    for (i in 1:2) if(i == 1){ plot(p_B) }else{plot(p_S)}
    devAskNewPage(options("device.ask.default")[[1]])
  }else{
    grid.arrange(p_B, p_S, ncol = 1)
  }
  
}