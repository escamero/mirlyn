#' Alpha Diversity Dataframe
#'
#' [alphadivDF()] will generate a dataframe of alpha-diversity values from the `mirl` object
#' allowing for characterization of uncertainty and variation introduced through random
#' subsampling during rarefying.
#'
#' @param x The `mirl` object.
#' @param diversity Diversity index to be applied. By default, the Shannon Index will
#'   be utilized in the generation of the dataframe. Diversity indexes available in vegan
#'   are supported.
#'
#' @return A `data.frame` consisting of sample metadata and selected diversity indexes.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#'
#' \dontrun{
#' mirlexample <- mirl(example, rep = 100)
#' alphadiv_df <- alphadivDF(mirlexample)
#' }
#'
#' @export
alphadivDF <- function(x, diversity="shannon"){
  md <- sample_data(x[[1]])
  div <- diversity(t(repotu_df(x)), index=diversity)
  final <- suppressWarnings(cbind(md, DiversityIndex = div))
  final
}

#' Visualization of Alpha Diversity Indexes
#'
#' [alphawichVis()] will visualize the diversity index values generated from [alphadivDF()].
#'
#' @param alphawichDF The `alphawichDF` object.
#' @param yvar The `yvar` column name. By default, the diversity index value will be plotted.
#' @param xvar The `xvar` column name. Users must specify the metadata column to be plotted
#'   on the x-axis.
#' @param colorvar The metadata column name that will be used for colour assignment.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' # This is one way of using it.
#' library(mirlyn)
#' data(example)
#'
#' \dontrun{
#' example <- mirl(example, rep = 100)
#' example <- alphadivDF(example)
#' alphawichVis(example, "Id")
#' }
#'
#' @export
alphawichVis <- function(alphawichDF, xvar = "Sample", yvar="DiversityIndex", colorvar=NULL) {
  if (is.null(colorvar)) {
    alpha <- ggplot(alphawichDF, aes_string(x = xvar, y = yvar))
  } else {
    alpha <- ggplot(alphawichDF, aes_string(x = xvar, y = yvar, color = colorvar))
  }
  alpha + theme_bw()+geom_jitter()+scale_y_continuous(limits=c(0, max(alphawichDF$DiversityIndex)*1.01))
}

