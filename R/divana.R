#' Diversity Analysis (DivAna)
#'
#' divana() will conduct a complete analysis from taxonomic composition bar charts to beta-diversity ordinations on samples. It includes the multiple iterations of rarefying and the generation of all associated plots with taxonomic composition, rarefaction curves, alpha diversity and beta-diversity analyses.
#'
#' @param x The phyloseq object.
#' @param xvar The independent categorical variable of interest.
#' @param sample The sample identifier.
#' @param lowabun The cutoff threshold for low abundance taxonomic clustering. By default, it will cluster all taxonomic groups present at less than 1% abundance into a low abundance group. Set to NULL if do not want any low abundance clustering.
#' @param reps The number of replicates to be included for rarefying. By default, it will conduct 1000 replicates.
#' @param alphadiv The alpha diversity metric. By default, the Shannon Index will be utilized in the generation of the dataframe. Diversity indexes available in vegan are supported.
#' @param taxcols A colour vector for taxonomic composition bar charts. Must include enough colours for all taxonomic groups present at all 7 levels.
#' @param divcols A colour vector for diversity analyses plots.
#' @param betatransform Supported standardization methods in community ecology in decostand() from library(vegan). Hellinger transformations must be applied to ecological data in order to allow for application in Euclidean space.
#' @param betadist The dissimilarity metric supported by vegdist() from library(vegan). By default, the Bray-Curtis dissimilarity matrix will be applied.
#'
#' @return A list of ggplot objects and a dataframe for alpha-diversity index values.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#'
#' \dontrun{
#' complete_analysis <- divana(example, "time")
#' }
#'
#' @export
divana <- function(x, xvar, sample = "Sample", lowabun = 0.01, reps = 1000, alphadiv = "Shannon", taxcols = NULL, divcols = NULL, betatransform = "hellinger", betadist = "bray"){
  #Taxonomic Bar Graph %
  fulltaxbar <- fullbartax(x, xvar, lowabun = lowabun, taxcols)
  #Rarefy Curve
  xrare <- rarefy_whole_rep(x, rep = reps)
  rarecurve <- rarecurve(xrare, sample = sample)
  #Multiple Rarefying
  multirare <- mirl(x)
  multirare_df <- repotu_df(multirare)
  #Alpha-wich
  alphawich_df <- alphadivDF(multirare_df)
  alphawich_plot <- alphawichVis(alphawich_df)
  #Alpha-Distributions
  alphacone_plot <- alphacone(x)
  #Beta-Diversity
  betamat <- betamatPCA(x, transformation = betatransform, dsim = betadist)
  betamat_vis <- betamatPCAvis(betamat, groups = xvar, reps = reps, colours = divcols)

  list(
    fulltaxbar = fulltaxbar,
    rarecurve = rarecurve,
    alphawich_df = alphawich_df,
    alphawich_plot = alphawich_plot,
    alphacone_plot = alphacone_plot,
    betamat_vis = betamat_vis
  )

}
