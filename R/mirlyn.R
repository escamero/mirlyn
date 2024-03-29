#' mirlyn: A package for library size normalization using repeated iterations of rarefaction
#'
#' @description
#'
#' `mirlyn` is a package that allows for library size normalization for
#' diversity analysis. This package also includes supplementary functions to
#' generate taxonomic composition bar graphs and the assignment of ASV names to
#' unique and easily referenced identities. The majority of the functions
#' included in mirlyn serve to account for the issue of disparity in liblrary
#' sizes between amplicon sequencing samples. The variation seen between library
#' sizes of different amplicon sequencing samples is not representative of true
#' biological variation and requires a normalization technique to allow for
#' sample comparison in diversity analysis without inherent bias. While a variety
#' of library normalization techniques have been proposed in the literature,
#' mirlyn repeatedly rarefies samples to characterize the sequence data lost in
#' the random subsampling process. The characterization of variation through
#' rarefying allows for the creation of comprehensive computational replicates
#' that provide stronger representation of the original biological communities
#' while also accounting for differences in library size.
#'
#' @docType package
#' @name mirlyn-pkg
#'
#' @import phyloseq
#' @import ggplot2
#' @importFrom factoextra fviz_pca_ind
#' @importFrom stats prcomp quantile p.adjust pnorm sd
#' @importFrom vegan vegdist decostand
#' @importFrom dplyr %>% group_by summarise filter distinct mutate_if funs
#' @importFrom parallel mclapply
#' @importFrom vegan diversity
#' @importFrom reshape2 melt
#' @importFrom tidyr gather pivot_wider
NULL
