divana <- function(x, sample = "SampleID", lowabun = 0.01, xvar, reps = 1000, alphadiv = "Shannon", taxcols, divcols, betatransform = "hellinger", betadist = "bray"){
  #Taxonomic Bar Graph %
  fulltaxbar <- fullbartax(x, lowabun = lowabun, xvar, cols)
  fulltaxbar
  #Rarefy Curve
  rarecurve <- rarefy_curve(x)
  rarecurve
  #Multiple Rarefying
  multirare <- mirl(x)
  multirare_df <- functionA(multirare)
  #Alpha-wich
  alphawich_df <- alphawichDF(multirare_df)
  alphawich_plot <- alphawichVis(alphawich_df)
  alphawich_pplot
  #Alpha-Distributions
  alphacone_plot <- alphacone(x)
  alphacone_plot
  #Beta-Diversity
  betamat <- betamatPCA(x, transformation = betatransform, distance = distance)
  betamat_vis <- betamatPCAvis(betamat, groups = xvar, reps = reps, cols = divcols)
  betamat_vis
}
