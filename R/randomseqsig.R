createobj <- function(x){
  obj <- otu_table(x)
  obj <- obj[rowSums(obj) > 0, ]
  obj
}

shuffle <- function(obj, randomx){
  randomx <- apply(obj, 2, sample)
  randomx
}

taxstruc <- function(x, obj,  taxlevel = "Phylum"){
  taxstruc <- structure(as.character(tax_table(x)[,taxlevel]), names = rownames(tax_table(x)))
  taxstruc <- taxstruc[rownames(obj)]
  taxstruc
}


shuf_and_count <- function(x, logivec) {
  x <- apply(x, 2, sample)
  colSums(x[logivec, ]) / colSums(x)
}


objshuffle <- function(obj, taxstruc, group ="Cyanobacteria", nshuff = 1000){
  objshuff <- lapply(1:nshuff, function(x) shuf_and_count(obj, taxstruc == group))
  objshuff <- do.call(rbind, objshuff)
  rownames(objshuff) <- 1:nshuff
  objshuff <- reshape2::melt(objshuff)
  objshuff
}

#real-mean/standard dev
#1/result squared

chebyshevstat <- function(stats){
  real_mean_stdev_calc <- (stats$real-stats$mean)/stats$sd
  overresultsqrd <- 1/(real_mean_stdev_calc^2)
  stats$pvalue <- overresultsqrd
  stats$adjpvalue <- p.adjust(overresultsqrd, method = "bonferroni" )
  stats
}

zscorestat <- function(stats) {
  res <- (stats$real - stats$mean) / stats$sd
  stats$pvalue <- pnorm(-abs(res))
  stats$adjpvalue <- p.adjust(stats$pvalue, method = "bonferroni")
  stats
}



get_stats <- function(x) {
  x %>%
    group_by(Var2) %>%
    summarise(mean = mean(value), sd = sd(value))
}

get_real_comp <- function(obj, logivec) {
  colSums(obj[logivec, ]) / colSums(obj)
}

zscorestat <- function(stats) {
  res <- (stats$real - stats$mean) / stats$sd
  stats$pvalue <- pnorm(-abs(res))
  stats$adjpvalue <- p.adjust(stats$pvalue, method = "bonferroni")
  stats
}

#' Compositional Significance Testing
#'
#' The `randomseqsig()` function will identify whether a taxonomic group of interest is significantly over-represented or underrepresented in the community. 
#'
#' @param x is a `phyloseq` object.
#' @param taxlevel is the taxonomic level of the group of interest.
#' @param group is the name of the taxonomic group of interest. 
#' @param nshuff is the number of times data shuffling is performed. 
#'
#' @return a list object which contatins: 1) `objshuff`: The shuffled sequence variant counts. 2) `stats`: The mean, standard deviation, real counts and 
#' p-values of the shuffled data.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#' significance_example <- randomseqsig(example, nshuff = 10)
#'
#' @export
randomseqsig <- function(x, taxlevel = "Phylum", group = "Cyanobacteria", nshuff = 1000){
  obj <- createobj(x)
  taxstruc <- taxstruc(x, obj, taxlevel)
  objshuff <- objshuffle(obj, taxstruc, group = group, nshuff = nshuff)
  obj_stats <- get_stats(objshuff)
  real_abun <- get_real_comp(obj, taxstruc == group)
  obj_stats$real <- real_abun
  obj_stats <- zscorestat(obj_stats)
  list(objshuff = objshuff, stats = obj_stats)
}
