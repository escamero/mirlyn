

rarefy_whole <- function(x, steps = seq(from = 0.001, to = 1, by = 0.01)){
  libsize <- sample_sums(x)
  libsizes <- max(libsize) * steps
  samplenames <- colnames(otu_table(x))
  #libsizerare <- seq(from = 0, to = libsize, by = 10)
  rarefy <- lapply(libsizes, function(y) rarefy_even_depth(x, y, verbose = FALSE))
  rarefy <- lapply(rarefy, function(y) apply(otu_table(y), 2, function(z) sum(z != 0)))
  rarefy <- lapply(rarefy, function(y) data.frame(Sample = samplenames, ObsASVCount = unname(y[samplenames]), row.names = NULL))
  rarefy <- mapply(function(y, z) { y$LibSize <- z; y }, rarefy, libsizes, SIMPLIFY = FALSE)
  rarefy <- as.data.frame(do.call(rbind, rarefy))
  rarefy
}

#' Rarefy Entire Sample Set
#'
#' rarefy_whole_rep() will rarefy each sample for different library sizes repeatedly in order to generate a representative set of data to be used in rarecurve() for generating rarefaction curves.
#'
#' @param x The phyloseq object.
#' @param rep The number of replicates to be performed. More replicates will take longer to conduct but will provide a smoother distribution for analyses.
#' @param steps The steps in library sizes to rarefy to.
#' @param set.seed The seed value for reproducibility.
#' @param mc.cores From [parallel::mclapply()].
#' @param intercept Create data points at (0,0).
#'
#' @return A dataframe containing the observed sequence counts for samples of varying rarefied library sizes.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#'
#' \dontrun{
#' rarefy_whole_rep_ex <- rarefy_whole_rep(example, rep = 100)
#' }
#'
#' @export
rarefy_whole_rep <- function(x, rep = 100, steps = seq(from = 0.001, to = 1, by = 0.01), set.seed = NULL, mc.cores = 1L, intercept = TRUE){
  if (!is.null(set.seed)) set.seed(set.seed)
  libsizes <- sample_sums(x)
  libsizes <- max(libsizes) * steps
  samplenames <- colnames(otu_table(x))
  rarefy_rep <- mclapply(seq_len(rep), function(y) suppressMessages(rarefy_whole(x, steps = steps)), mc.cores = mc.cores)
  all_Obs <- lapply(rarefy_rep, function(y) as.matrix(y[, 2, drop = FALSE]))
  all_Obs <- do.call(cbind, all_Obs)
  all_Obs <- rowMeans(all_Obs)
  rarefy_rep <- rarefy_rep[[1]]
  rarefy_rep$ObsASVCount <- all_Obs
  if (intercept) {
    zeroLibSize <- data.frame(Sample = unique(rarefy_rep$Sample), ObsASVCount = 0, LibSize = 0, row.names = NULL)
    rbind(rarefy_rep, zeroLibSize)
  } else {
    rarefy_rep
  }
}

#' Rarefaction Curve
#'
#' rarecurve() will generate a visualization of the rarefaction curve showing the associated library size.
#'
#' @param x The rarefy_whole_rep() object.
#' @param sample The sample ID which will be used to assign colours to different samples. Leave blank if don't want colours.
#'
#' @return A ggplot object of the rarefaction curve showing the observed ASV count according to library size of each sample.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#'
#' \dontrun{
#' rarefy_whole_rep_ex <- rarefy_whole_rep(example, rep = 100)
#' rarecuruve_ex <- rarecurve(rarefy_whole_rep_ex, sample = "Sample")
#' }
#'
#' @export
rarecurve <- function(x, sample = "Sample"){
  if (is.null(sample)){
    curve <- ggplot(x, aes_string(x = "LibSize", y = "ObsASVCount"))
    }else{
    curve <- ggplot(x, aes_string(x = "LibSize", y = "ObsASVCount", colour = sample))
  }
  curve + geom_line()+theme_bw()+scale_y_continuous(limits=c(0, max(x$ObsASVCount)*1.01))
}


