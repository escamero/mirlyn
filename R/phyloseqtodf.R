getmeta <- function(x){
  metadata <- as.matrix(x@sam_data)
  sample <- rownames(metadata)
  metadata <- cbind(sample, metadata)
  metadata
}

getotu <- function(x){
  otutable <- as.data.frame(x@otu_table@.Data)
  otutable
}

gettax <- function(x){
  taxtable <- as.data.frame(x@tax_table@.Data)
  taxtable
}

getotutax <- function(x){
  otutaxtable <- merge(getotu(x), gettax(x), by = "row.names", all.x=TRUE)
  otutaxtable <- gather(otutaxtable, "sample", "abundance", 1 + 1:ncol(getotu(x)))
  otutaxtable
}

#' Conversion of `phyloseq` Object to Dataframe
#'
#' `phyloseq_to_df()` will create a dataframe including OTU abundances, taxonomic classification and metadata extracted from the `phyloseq` object. 
#' This dataframe can be exported as a CSV file for easy viewing. 
#'
#'
#' @param x The `phyloseq` object.
#'
#' @return A dataframe object.
#' 
#' #' @examples
#' library(mirlyn)
#' data(example)
#' example_df <- phyloseq_to_df(example)
#' 
#' @export
phyloseq_to_df <- function(x){
  physeqdf <- merge(getotutax(x), getmeta(x))
  physeqdf
}
