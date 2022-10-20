#' Global
#' 
#' Declare global variables
#' @name global
#' @importFrom utils globalVariables
NULL

if(getRversion() >= "2.15.1")  utils::globalVariables(c("x", "Predmean", "Pred025", "Pred975"), add = F)