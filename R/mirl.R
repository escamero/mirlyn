# Used to calculate library sizes
sample_sums <- phyloseq::sample_sums
#' Multiple Iterations of Rarefying Libraries
#' 
#' mirl() will repeatedly rarefy to a user defined library size, for a specified number of replications to characterize the variation introduced through random subsampling.
#'
#' @param x The phyloseq object
#' @param libsize The specified library size to rarefy to. By default, mirl() will rarefy to the minimum library size found in samples. 
#' @param rep The number of times rarefying will be repeated. By default, mirl() will repeatedly rarefy 1000 times. 
#' @param set.seed The seed value for reproducibility. 
#' @param trimOTUs A Boolean value determining whether OTUs that are absent from all samples after rarefying should be removed. 
#' @param replace A Boolean value determining whether subsampling should be performed with or without replacement. See Cameron et al., 2020 for details. 
#' @examples
#' library(mirlyn)
#' mirl_object <- mirl(x, libsize = 10000, rep = 1000, set.seed=120)
#' @export
mirl <- function(x, libsize=min(sample_sums(x)), rep=1000, set.seed=NULL, trimOTUs=FALSE, replace=FALSE){
  if (!is.null(set.seed)) set.seed(set.seed)
  replicate(rep, phyloseq::rarefy_even_depth(x, samplesize=libsize, trimOTUs = trimOTUs, replace = replace))
}
