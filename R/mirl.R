# Used to calculate library sizes
# sample_sums <- phyloseq::sample_sums

#' Multiple Iterations of Rarefying Libraries
#'
#' [mirl()] will repeatedly rarefy to a user defined library size, for a specified number
#' of replications to characterize the variation introduced through random subsampling.
#' The implementation of mirl() allows for rarefying to be a statistically valid library
#' size normalization technique for diversity analyses. Users are encouraged to conduct
#' exploratory analysis to identify optimal rarefied library sizes for their data that
#' minimizes variation in samples. Users may have to choose between rarefying to smaller
#' than preferred library sizes or discarding small library size samples but the
#' implementation of [mirl()] generally allows for rarefying to small library sizes at the
#' loss of resolution and precision in results.
#'
#' @param x The `phyloseq` object
#' @param libsize The specified library size to rarefy to. By default, [mirl()] will rarefy
#'   to the minimum library size found in samples.
#' @param rep The number of times rarefying will be repeated. By default, [mirl()] will
#'   repeatedly rarefy 1000 times.
#' @param set.seed The seed value for reproducibility.
#' @param trimOTUs A Boolean value determining whether OTUs that are absent from all
#'   samples after rarefying should be removed.
#' @param replace A Boolean value determining whether subsampling should be performed
#'   with or without replacement. See Cameron et al., 2020 for further details.
#' @param mc.cores From [parallel::mclapply()].
#'
#' @return A list of `phyloseq` objects.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#'
#' \dontrun{
#' mirl_object <- mirl(example, libsize = 10000, rep = 100, set.seed=120)
#' }
#'
#' @export
mirl <- function(x, libsize=min(phyloseq::sample_sums(x)), rep=1000, set.seed=NULL, trimOTUs=FALSE, replace=FALSE, mc.cores = 1L){
  if (!phyloseq::taxa_are_rows(x))
    stop("The mirlyn package requires phylsoseq::taxa_are_rows(x) to be TRUE.")
  if (!is.null(set.seed)) set.seed(set.seed)
  set.seed <- sample.int(1e9, rep)
  mirlobj <- mclapply(seq_len(rep),
    function(y) suppressMessages(rarefy_even_depth(x, sample.size=libsize, trimOTUs = trimOTUs,
        replace = replace, verbose = FALSE, rngseed = set.seed[y])),
    mc.cores = mc.cores)
  l1 <- length(sample_names(x))
  l2 <- length(sample_names(mirlobj[[1]]))
  if (l1 > l2) {
    warning(l1 - l2, " samples were dropped due to too ",
      "few counts.")
  }
  mirlobj
}
