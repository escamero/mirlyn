#' Heatmaps for Taxonomic Group of Interest
#'
#' The `plot_heat()` function will create heatmaps visualizing the relative abundance of a taxonomic group of interest from a dataframe.
#'
#' @param x is a `phyloseqtodf` object.
#' @param taxlevel is the taxonomic level of the taxa you're interested in. By default, it will look at the phylum level.
#' @param taxaname is the taxonomic name of interest.
#' @param xvar is a categorical variable selected from the metadata.
#' @param yvar is a secondary categorical variable selected from the metadata.
#' @param fillvar is the column name of the abundance.
#'
#' @return a `ggplot` object.
#'
#' @examples
#' library(mirlyn)
#' library(dplyr)
#' data(example)
#'
#' example_df <- phyloseq_to_df(example)
#' example_df_phylum <- example_df %>%
#'     group_by(sample, Id, Phylum) %>%
#'     summarise(abundance = sum(abundance))%>%
#'     mutate(Proportion = abundance/sum(abundance)*100)
#' plot_heat(example_df_phylum, taxlevel = "Phylum",taxaname = "Cyanobacteria",
#'     xvar = "sample", yvar = "Id", fillvar = "Proportion")
#'
#' @export
plot_heat <- function(x, taxlevel,taxaname, xvar, yvar, fillvar){
  tax_x <- x %>% filter(taxlevel == taxaname)
  heatmap <-ggplot(x, aes_string(xvar, yvar, fill= fillvar))+geom_tile()+theme_bw()+scale_y_discrete(expand = c(0,0))+scale_x_discrete(expand = c (0,0))
  heatmap
}


