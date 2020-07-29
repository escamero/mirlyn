functionRep <- function(x, libs = c(15000, 16000, 17000, 18000, 19000), rep = 1000, set.seed = NULL) {
  if (!is.null(set.seed)) set.seed(set.seed)
  out <- vector("list", length(libs))
  for (i in seq_along(out)) {
    message("Performing replicates for library size of ", libs[i], "...")
    out[[i]] <- replicate(rep, phyloseq::rarefy_even_depth(x, libs[i], verbose = FALSE))
  }
  names(out) <- libs
  out
}

functionA <- function(x) {
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

functionA2 <- function(x) {
  lapply(x, functionA)
}

functionB <- function(x, diversity = "shannon", quantiles = c(0.025, 0.975)) {
  x2 <- vegan::diversity(t(x), index = diversity)
  quantile(x2, quantiles)
}
functionB2 <- function(x, diversity = "shannon", quantiles = c(0.025, 0.975)) {
  y <- lapply(x, function(z) functionB(z, diversity, quantiles))
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
#' alphacone() will generate a distribution plot of the alpha diversity index for different library sizes. Different rarefied library sizes may have impact on both the measured value of the diversity index and the perceived variation in the index values. Generations of the distribution of the alpha diversity index allows for a comprehensive examination of the alpha diversity index values to account for variation introduced into the diversity index as an artifact of rarefied libraray size.
#'
#'

#' @export
alphacone <- function(obj, libs = c(15000, 16000, 17000, 18000, 19000), rep = 1000, set.seed = NULL, diversity = "shannon", lower.q = 0.025, upper.q = 0.975) {
  the_reps <- functionRep(obj, libs = libs, rep = rep, set.seed = set.seed)
  the_reps_fixed <- functionA2(the_reps)
  the_reps_q <- functionB2(the_reps_fixed, diversity = diversity, quantiles = c(lower.q, upper.q))
  the_reps_q
}
