#Lowabundance_tax clumps low abundance taxa into a "Low abundance" group. Default is set to 1%. 
lowabundance_tax <- function(x, taxrank = "Genus", lower_limit = 0.01) {
  x <- microbiome::transform(x, "compositional")
  x <- phyloseq::tax_glom(x, taxrank = taxrank)
  x <- phyloseq::psmelt(x)
  x[[taxrank]] <- as.character(x[[taxrank]])
  x[[taxrank]][x[["Abundance"]] < lower_limit] <- "<1% abundance"
  x
}

#If not using low abundance, will just use the raw sequence count values extracted from the phyloseq object using psmelt. 
raw_tax <- function(x, taxrank = "Genus"){
  x <- phyloseq::tax_glom(x, taxrank = taxrank)
  x <- phyloseq::psmelt(x)
  x
}

#' Taxonomic Composition Stacked Barcharts
#'
#' The bartax() function will create taxonomic composition barcharts. If desired, low-abundant taxonomic groups will be clustered into a single low abundant group. 
#' 
#' @param x is a phyloseq object.
#' @param lowabun is the cut off value for a low-abundant group. 
#' @param xvar is a categorical variable selected from the metadata.
#' @param yvar is the relative abundance of reads. 
#' @param taxrank is the specified taxonomic rank. 
#' @param cols is a vector containing specified colour values (sufficient number for number of taxonomic groups). 
#'
#' @return a ggplot object.
#' 
#' @examples 
#' library(mirlyn)
#' data(example)
#' cols <- c("black","darkgoldenrod1","dodgerblue","deeppink4","chartreuse3","burlywood4","navy","blueviolet", "tan2","lavenderblush3")
#' bartax(example, lowabun = 0.01, SampleID, yvar = "Abundance", taxrank = "Phylum", cols)
#'
#' @export
bartax <- function(x, lowabun = 0.01, xvar, yvar = "Abundance", taxrank = "Genus", cols){
  if (is.null(lowabun)){
    x <- lowabundance_tax(x, taxrank = taxrank, lower_limit = lowabun)
  }else{
    x <- raw_tax(x, taxrank = taxrank)
  }
  if(is.null(cols)){
    taxbar<- ggplot(x, aes(x=xvar, y=yvar, fill=taxrank))+geom_bar(stat="identity")+theme_bw()+scale_y_continuous(labels=function(x) paste0(x*100,"%"),limits=c(0,1.05),expand=c(0,0))
  }else{
  taxbar<- ggplot(x, aes(x=xvar, y=yvar, fill=taxrank))+geom_bar(stat="identity")+theme_bw()+scale_fill_manual(values=cols)+scale_y_continuous(labels=function(x) paste0(x*100,"%"),limits=c(0,1.05),expand=c(0,0))}
  taxbar
}
#' Taxonomic Composition Stacked Barcharts - 7 levels
#' 
#' The fullbartax() function will create taxonomic composition for barcharts for all 7 taxonomic levels.
#' @param x is a phyloseq object.
#' @param lowabun is the cut off value for a low-abundant group. 
#' @param xvar is a categorical variable selected from the metadata.
#' @param yvar is the relative abundance of reads. 
#' @param cols is a vector containing specified colour values (sufficient number for number of taxonomic groups found at all levels).
#' 
#' @examples
#' library(mirlyn)
#' data(example)   
#' fullbartax(example, lowabun = 0.01, xvar, yvar="Abundance")
#'                                                                                                                                               
#' @export
fullbartax <- function(x, lowabun = 0.01, xvar, yvar = "Abundance", cols){
  kingdom_plot <- bartax(x, xvar = xvar, taxrank = "Kingdom", cols = cols)
  phylum_plot <- bartax(x, xvar = xvar, taxrank = "Phylum", cols = cols)
  class_plot <- bartax(x,xvar = xvar, taxrank = "Class", cols = cols)
  order_plot <- bartax(x, xvar = xvar, taxrank = "Order", cols = cols)
  family_plot <- bartax(x, xvar = xvar, taxrank = "Family", cols = cols)
  genus_plot <- bartax(x, xvar = xvar, taxrank = "Genus", cols = cols)
  species_plot <- bartax(x, xvar = xvar, taxrank = "Species", cols = cols)
  kingdom_plot
  phylum_plot
  class_plot
  order_plot
  family_plot
  genus_plot
  species_plot
}





