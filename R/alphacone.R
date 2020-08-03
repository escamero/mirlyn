

rarefy_rep_otu <- function(x, rep = 100, steps = seq(from = 0.001, to = 1, by = 0.01), set.seed = NULL){
  if (!is.null(set.seed)) set.seed(set.seed)
  rarefy_rep_otu <- lapply(seq_len(rep), function(y) rarefy_whole(x, steps = steps))
  rarefy_rep_otu
}
#Old function for alphacone
#functionRep <- function(x, libs = c(15000, 16000, 17000, 18000, 19000), rep = 1000, set.seed = NULL) {
 # if (!is.null(set.seed)) set.seed(set.seed)
  #out <- vector("list", length(libs))
  #for (i in seq_along(out)) {
   # message("Performing replicates for library size of ", libs[i], "...")
    #out[[i]] <- replicate(rep, rarefy_even_depth(x, libs[i], verbose = FALSE))
  #}
  #names(out) <- libs
  #out
#}

#Replicate OTU Table to Dataframe - Single Library Size
repotu_df <- function(x) {
  newdat <- lapply(x, otu_table)
  for (i in seq_along(newdat)) {
    colnames(newdat[[i]]) <- paste0(i, "-", colnames(newdat[[i]]))
  }
  newdat2 <- lapply(newdat, as.data.frame)
  for (i in seq_along(newdat2)) {
    newdat2[[i]]$SeqId <- rownames(newdat2[[i]])
    rownames(newdat2[[i]]) <- NULL
  }
  finaldat <- newdat2[[1]]
  for (i in seq_along(newdat2)[-1]) {
    finaldat <- merge(finaldat, newdat2[[i]], by = "SeqId", all = FALSE)
  }
  rownames(finaldat) <- finaldat$SeqId
  finaldat[, -1]
}

#Replicate OTU Table to Dataframe - Multiple Library Sizes
repotu_libsize_df <- function(x) {
  lapply(x, repotu_df)
}

div_quantile <- function(x, diversity = "shannon", quantiles = c(0.025, 0.975)) {
  x2 <- vegan::diversity(t(x), index = diversity)
  quantile(x2, quantiles)
}
div_quantile_df <- function(x, diversity = "shannon", quantiles = c(0.025, 0.975)) {
  y <- lapply(x, function(z) repotu_libsize_df(z, diversity, quantiles))
  for (i in seq_along(y)) {
    y[[i]] <- data.frame(Lower = y[[i]][1], Upper = y[[i]][2])
  }
  y2 <- do.call(rbind, y)
  y2$LibSize <- as.numeric(names(x))
  rownames(y2) <- NULL
  y2
}

#' Alpha Diversity Index Distributions
#'
#' alphacone() will generate a distribution plot of the alpha diversity index for different library sizes. Different rarefied library sizes may have impact on both the measured value of the diversity index and the perceived variation in the index values. Generations of the distribution of the alpha diversity index allows for a comprehensive examination of the alpha diversity index values to account for variation introduced into the diversity index as an artifact of rarefied library size.
#' @param x The phyloseq object.
#' @param rep The number of replications of rarefying to perform.
#' @param steps The library sizes to rarefy to as a function of the library sizes of samples.
#' @param diversity The diversity index. By default, the Shannon Index will be utilized in the generation of the dataframe. Diversity indexes available in vegan are supported.
#' @param lower.q The lower quantile for the distribution cone. By default it will provide the 2.5th.
#' @param upper.q The upper quantile for the distribution cone. By default, the 97.5th will be applied.
#'
#' @return A ggplot object showing the distribution of the diversity index.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#'
#' alphacone_example <- alphacone(example, rep = 10)
#'
#'
#' # Don't think this works on multiple samples right now.
#' @export
alphacone <- function(x, rep = 1000, steps = seq(from = 0.001, to = 1, by = 0.01), diversity = "shannon", lower.q = 0.025, upper.q = 0.975) {
  the_reps <- rarefy_rep_otu(x, rep = rep, steps = steps)
  the_reps_fixed <- repotu_libsize_df(the_reps)
  the_reps_q <- div_quantile_df(the_reps_fixed, diversity = diversity, quantiles = c(lower.q, upper.q))
  the_reps_q
}
