rarefy_rep_otu <- function(x, rep = 100, steps = seq(from = 0.001, to = 1, by = 0.01), set.seed = NULL){
  if (!is.null(set.seed)) set.seed(set.seed)
  rarefy_rep_otu <- lapply(seq_len(rep), function(y) rarefy_whole(x, steps = steps))
  rarefy_rep_otu
}

rarefy_lib_otu <- function(x, libsizes) {
  lapply(libsizes, function(y) otu_table(rarefy_even_depth(x, y, verbose = FALSE)))
}

rarefy_lib_otu_rep <- function(x, libsizes, rep) {
  lapply(seq_len(rep), function(y) rarefy_lib_otu(x, libsizes))
}

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
  y <- lapply(x, function(z) div_quantile(z, diversity, quantiles))
  for (i in seq_along(y)) {
    y[[i]] <- data.frame(Lower = y[[i]][1], Upper = y[[i]][2])
  }
  y2 <- do.call(rbind, y)
  y2$LibSize <- as.numeric(names(x))
  rownames(y2) <- NULL
  y2
}

#' Alpha Diversity Metric Distributions
#'
#' [alphacone()] will generate a dataframe of the alpha diversity metric for different
#' library sizes. Different rarefied library sizes may have impact on both the measured
#' value of the diversity metric and the perceived variation in the metric values.
#' Generations of the distribution of the alpha diversity metric allows for a comprehensive
#' examination of the alpha diversity metric values to account for variation introduced into
#' the diversity metric as an artifact of rarefied library size.
#'
#' @param x The `phyloseq` object.
#' @param rep The number of replications of rarefying to perform.
#' @param steps The library sizes to rarefy to as a function of the library sizes of samples.
#' @param diversity The diversity metric. By default, the Shannon Index will be utilized in
#'   the generation of the dataframe. Diversity indexes available in vegan are supported.
#' @param lower.q The lower quantile for the distribution cone. By default it will provide the 2.5th.
#' @param upper.q The upper quantile for the distribution cone. By default, the 97.5th will be applied.
#' @param replace Whether to use replacement during sampling.
#'
#' @return A `data.frame` of the distribution of the diversity metric.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#'
#' \dontrun{
#' alphacone_example <- alphacone(example, rep = 100)
#' }
#'
#'
#' @export
alphacone <- function(x, rep = 1000, steps = seq(from = 0.001, to = 1, by = 0.01),
  diversity = "shannon", lower.q = 0.025, upper.q = 0.975, replace = FALSE) {
  libsizes <- max(sample_sums(x)) * steps
  samplenames <- sample_names(x)
  meta <- sample_data(x)
  bylibsize <- vector("list", length(libsizes))
  for (i in seq_along(bylibsize)) {
    bylibsize[[i]] <- lapply(seq_len(rep),
      function(y) suppressMessages(rarefy_even_depth(x, libsizes[i], verbose = FALSE, replace = replace)))
  }
  for (i in seq_along(bylibsize)) {
    for (j in seq_along(bylibsize[[i]])) {
      otutable <- t(otu_table(bylibsize[[i]][[j]]))
      if (nrow(otutable) == 1) {
        otutable <- as.data.frame(otutable)
        tmpnames <- rownames(otutable)
      }
      bylibsize[[i]][[j]] <- diversity(otutable, index = diversity)
      if (nrow(otutable) == 1) {
        names(bylibsize[[i]][[j]]) <- tmpnames
      }
    }
  }
  for (i in seq_along(bylibsize)) {
    for (j in seq_along(bylibsize[[i]])) {
      bylibsize[[i]][[j]] <- data.frame(
        Sample = names(bylibsize[[i]][[j]]), LibSize = libsizes[i],
        Rep = j, DiversityIndex = bylibsize[[i]][[j]])
    }
    bylibsize[[i]] <- do.call(rbind, bylibsize[[i]])
  }
  bylibsize <- do.call(rbind, bylibsize)
  rownames(bylibsize) <- NULL
  bylibsize <- bylibsize %>% group_by(Sample, LibSize) %>%
    summarise(LowerQ = quantile(DiversityIndex, lower.q),
      UpperQ = quantile(DiversityIndex, upper.q))
  bylibsize <- as.data.frame(bylibsize)
  bylibsize <- cbind(bylibsize, meta[bylibsize$Sample, ])
  rownames(bylibsize) <- NULL
  bylibsize
}

#' Alpha Diversity Metric Distributions Visualization
#'
#' [alphaconeVis()] will generate a visualization of the the distributions of the alpha
#' diversity metric `data.frame` generated using [alphacone()].
#'
#' @param x The [alphacone()] dataframe
#' @param cols The metadata column that will be used for colour assignment
#' @param alpha The opacity of the geom_ribbon. By default, a value of 0.7 will be applied.
#'   For full opacity, set alpha = 1, for more transparency, decrease alpha accordingly.
#'
#' @return a `ggplot` object.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#'
#' \dontrun{
#' alphacone_example <- alphacone(example, rep = 100)
#' alphaconeVis(alphacone_example, cols = "Sample")
#' }
#'
#' @export
alphaconeVis <- function(x, cols = "Sample", alpha = 0.7){
  alphaconevis <- ggplot(x, aes(LibSize)) +
    geom_ribbon(aes_string(ymin = "LowerQ", ymax = "UpperQ", fill = cols),
      alpha = alpha) +
    theme_bw()
  alphaconevis
}





