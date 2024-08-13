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
alphadivDF <- function(x, diversity = c("shannon", "simpson", "invsimpson")) {
    md <- sample_data(x[[1]])
    
    # Calculate observed richness
    t_otu_table <- t(rep_otu_df(x))
    observed <- specnumber(t_otu_table)
    observed_richness_df <- rownames_to_column(as.data.frame(observed), var = "Unique_ID")
    
    div_df_list <- lapply(diversity, function(index) {
        div_values <- vegan::diversity(t(rep_otu_df(x)), index = index)
        if (!is.matrix(div_values)) {
            div_values <- as.matrix(div_values, ncol = 1, dimnames = list(rownames(div_values), NULL))
        }
        colnames(div_values) <- paste0(index, "_", colnames(div_values))
        return(data.frame(Unique_ID = rownames(div_values), div_values, stringsAsFactors = FALSE))
    })
    
    final <- reduce(div_df_list, function(df1, df2) inner_join(df1, df2, by = "Unique_ID"))
    
    # Merge observed richness with the final dataframe
    final <- merge(final, observed_richness_df, by = "Unique_ID", all.x = TRUE)
    
    final <- cbind(md, final)
    
    return(final)
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

