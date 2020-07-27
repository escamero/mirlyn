#' @export
alphawichDF <- function(x, diversity="shannon"){
  md <- sam_data(x[[1]])
  div <- vegan::diversity(t(functionA(x)), index=diversity)
  final <- cbind(md, DiversityIndex = div)
  final
}

#' This is the alphawhichvis function.
#'
#' The description.
#'
#' @param alphawichDF The alphawichDF object.
#' @param yvar The yvar column name.
#' @param xvar The xvar column name.
#'
#' @return A ggplot object.
#'
#' @examples
#' # This is one way of using it.
#' #library(mirlyn)
#' #data(exampleobject)
#'
#' #alphawichvis(exampleobject)
#'
#' @export
alphawichVis <- function(alphawichDF, yvar="DiversityIndex", xvar, colorvar="Group") {
  ggplot(alphawichDF, aes_string(x = xvar, y = yvar, color=colorvar))+theme_bw()+geom_jitter()+scale_y_continuous(limits=c(0, max(alphawichDF$DiversityIndex)*1.01))
}

