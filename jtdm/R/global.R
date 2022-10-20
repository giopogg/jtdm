#' Global
#' 
#' Declare global variables
#' @importFrom utils globalVariables
NULL

## quiets concerns of R CMD check re: the .'s that appear in pipelines
if(getRversion() >= "2.15.1")  utils::globalVariables(c("x", "Predmean", "Pred025", "Pred975"))