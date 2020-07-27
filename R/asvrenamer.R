#' Assign ASV identifiers to sequence variants
#'
#' The asv_rename function will assign unique ASV identifiers that are easy to each sequence variant and will write a tsv file to your working directory for a reference list of all ASV names.
#' @param x The phyloseq object.
#' @param path The path to the tsv generated file of sequence variant identifiers.
#' 
#' @return x the modified phyloseq object.
#' @export
asv_rename <- function(x, path = getwd()){
  treefix <- x@phy_tree$tip.labels
  fixed <- paste0("ASV", seq_along(tofix))
  names(fixed) <- tofix
  asvdf <- data.frame(variant = names(fixed), newname = fixed, stringsAsFactors = FALSE)
  readr::write_tsv(asvdf, path/ASV_reflist.tsv)
}

#' Assign ASV identifiers to FASTA
#'
#' The filt_fasta_rename function will assign corresponding ASV identifiers to a FASTA file for use in other applications. 
#' @param filename The FASTA file to assign ASV identifiers to.
#' @param fixed The ASV identifiers.
#' 
#' @return FASTA file with modified ASV identifiers.
#' @export
filt_fasta_rename <- function(filename, fixed) {
  d <- Biostrings::readDNAStringSet(filename)
  names(d)[names(d) %in% names(fixed)] <- fixed[names(d)[names(d) %in% names(fixed)]]
  Biostrings::writeXStringSet(d, filename)}
