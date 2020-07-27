

rarefy_whole <- function(x, step = 10){
  libsize <- sample_sums(x)
  libsizerare <- c(seq_along(from = 0, to = libsize, by = 10))
  rarefy <- rarefy_even_depth(x, libsizerare[i])
  rarefy
}

#' @export
rarecurve <- function(x, sample = "SampleID"){
  rarefywhole <- rarefy_whole(x)
  rarefywhole_df <- functionA(rarefywhole)
  variantnum <- nrows(rarefywhole_df)
  curve <- ggplot(rarefywhole_df, aes(x = libsize, y = variantnum, color = sample))
  curve
}


