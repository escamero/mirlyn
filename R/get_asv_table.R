strip_meta <- function(x, OTU = "OTU", Sample = "sample", Abundance = "abundance"){
  x <- x[c(OTU, Sample, Abundance)]
  x <- as.data.frame(pivot_wider(x, names_from = all_of(Sample), values_from = all_of(Abundance)))
  x
}

relbuns <- function(x){
  x <- x %>% mutate_if(is.numeric, funs((./sum(.))*100))
  x
}


get_taxid <- function(x, OTU = "OTU", Kingdom = "Kingdom", Phylum = "Classified.Phylum", Class = "Classified.Class", Order = "Classified.Order", Family = "Classified.Family", Genus = "Classified.Genus", Species = "Classified.Species"){
  x <- x[c(OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species)]
  x <- distinct(x)
  x
}


#' Compiled ASV Table with Abundance & Taxonomic Classification
#'
#' `get_asv_table()` will produce a compiled ASV table with abundance and taxonomic classifications using the output from `phyloseqtodf()`.
#' In combination with the `filter()` function from `dplyr`, users may select a specific taxonomic group.
#'
#'
#' @param x The `phyloseqtodf` object
#' @param OTU The ASV column name. By default, the OTU column name is "OTU".
#' @param Sample The Sample ID column name. By default, the Sample ID column  name is "sample" from the `phyloseqtodf()` function.
#' @param Abundance The Abundance column name. By default, the Sample ID column  name is "abundance" from the `phyloseqtodf()` function.
#' @param Kingdom The taxonomic Kingdom column name. By default, the taxonomic Kingdom column name is "Kingdom" from the `phyloseqtodf()` function.
#' @param Phylum The taxonomic Phylum column name. By default, the taxonomic Phylum column name is "Classified.Phylum" from the `phyloseqtodf()` function.
#' @param Class The taxonomic Class column name. By default, the taxonomic Class column name is "Classified.Class" from the `phyloseqtodf()` function.
#' @param Order The taxonomic Order column name. By default, the taxonomic Order column name is "Classified.Order" from the `phyloseqtodf()` function.
#' @param Family The taxonomic Family column name. By default, the taxonomic Family column name is "Classified.Family" from the `phyloseqtodf()` function.
#' @param Genus The taxonomic Genus column name. By default, the taxonomic Genus column name is "Genus" from the `phyloseqtodf()` function.
#' @param Species The taxonomic Species column name. By default, the taxonomic Species column name is "Classified.Species" from the `phyloseqtodf()` function.
#' @param proportion A Boolean value determining whether abundances of ASV should be expressed as raw read counts or % abundance of the library. By default,
#' the % abundance will be applied.
#'
#' @return A dataframe.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#' example_df <- phyloseq_to_df(example)
#'
#' # ASV Table for All ASV
#' example_asv_table <- get_asv_table(example_df, OTU = "Row.names")
#'
#' # Filter by Taxonomic Level of Interest
#' library(dplyr)
#' cyano_example <- example_df %>% filter(Phylum == "Cyanobacteria")
#' cyano_example_asv_table <- get_asv_table(cyano_example, OTU = "Row.names")
#'
#' @export
get_asv_table <- function(x, OTU = "OTU", Sample = "sample", Abundance = "abundance", Kingdom = "Kingdom", Phylum = "Phylum", Class = "Class", Order = "Order", Family = "Family", Genus = "Genus", Species = "Species", proportion = TRUE){
  nometa_x <- strip_meta(x, OTU, Sample, Abundance)
  taxid_x <- get_taxid(x, OTU, Kingdom, Phylum, Class, Order, Family, Genus, Species)
  if(proportion){
    nometa_x <- relbuns(nometa_x)
  }
  merge(nometa_x, taxid_x)
}
