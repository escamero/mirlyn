#' Assign ASV identifiers to sequence variants
#'
#' The asv_rename function will assign unique ASV identifiers that are easy to each sequence variant and will write a tsv file to your working directory for a reference list of all ASV names.
#' @param x The phyloseq object.
#' @param file The path to the tsv generated file of sequence variant identifiers.
#'
#' @return x the modified phyloseq object.
#'
#' @examples
#' library(mirlyn)
#' library(phyloseq)
#' data(GlobalPatterns)
#' asv_rename(GlobalPatterns)
#'
#' @export
asv_rename <- function(x, file = file.path(getwd(), "ASVids.tsv")){
  treefix <- phy_tree(x)$tip.label
  fixed <- paste0("ASV", seq_along(treefix))
  names(fixed) <- treefix
  asvdf <- data.frame(Variant = names(fixed), ASVid = fixed,
    stringsAsFactors = FALSE, row.names = NULL)
  readr::write_tsv(asvdf, file)
  taxa_names(x) <- fixed[taxa_names(x)]
  x
}

#' Assign ASV identifiers to FASTA
#'
#' The filt_fasta_rename function will assign corresponding ASV identifiers to a FASTA file for use in other applications.
#' @param fasta The FASTA file to assign ASV identifiers to.
#' @param fixed The ASV reference list created by [asv_rename()].
#'
#' @return `NULL`. The FASTA file is modified in place.
#'
#' @export
filt_fasta_rename <- function(fasta, fixed = file.path(getwd(), "ASVids.tsv")) {
  d <- Biostrings::readDNAStringSet(fasta)
  fixed <- readr::read_tsv(fixed)
  fixed <- structure(fixed$ASVid, names = fixed$Variant)
  names(d)[names(d) %in% names(fixed)] <- fixed[names(d)[names(d) %in% names(fixed)]]
  Biostrings::writeXStringSet(d, fasta)
  invisible(NULL)
}
