#' Beta Diversity Matrix - PCA
#'
#'  betamatPCA() will perform a PCA on the mirl() object.
#'
#' @param x The mirl() object.
#' @param transformation Supported standardization methods in community ecology in decostand() from library(vegan). Hellinger transformations must be applied to ecological data in order to allow for application in Euclidean space.
#' @param dsim The dissimilarity metric supported by vegdist() from library(vegan). By default, the Bray-Curtis dissimilarity matrix will be applied.
#'
#' @return A prcomp object
#'
#' @examples
#' library(mirlyn)
#' data(example)
#' \dontrun{
#' mirlexample <- mirl(example, rep = 100)
#' betamatPCA_object <- betamatPCA(mirlexample, dsim = "jaccard")
#' }
#'
#' @export
betamatPCA <- function(x, transformation="hellinger", dsim="bray"){
  if (is.null(transformation)){
    dist <- vegdist(as.matrix(t(repotu_df(x))))
  } else {
    transform <- decostand(as.matrix(t(repotu_df(x))), transformation)
    dist <- vegdist(transform, method=dsim)
  }
  pca <- prcomp(dist)
  pca
}

#' Beta Diversity Visualization
#'
#' betamatPCAvis() will produce a visualization from the prcomp object generated from betamatPCA_object.
#'
#' @param x The betamatPCAvis object.
#' @param geom The geom style to be aplied to the plot.
#' @param groups A vector containing the metadata variable being explored.
#' @param reps The number of replicates conducted in the initial mirl().
#' @param colours A vector containing colour values.
#'
#' @return A ggplot object.
#'
#' @examples
#' library(mirlyn)
#' data(example)
#' mirlexample <- mirl(example, rep = 10)
#' betamatPCA_object <- betamatPCA(mirlexample)
#'
#' \dontrun{
#' betamatPCAvis(betamatPCA_object, groups = c("A", "B", "C", "D","E","F"), reps = 10,
#'    colours = c("#000000", "#E69F00", "#0072B2", "#009E73", "#F0E442", "#D55E00"))
#' }
#'
#' @export
betamatPCAvis <- function(x, geom="point", groups, reps, colours){
  pcaplot <- fviz_pca_ind(x, geom=geom, habillage=rep(groups, reps),palette=colours, invisible="quali")+theme_bw()
  pcaplot
}

