#' @export
betamatPCA <- function(x, transformation="hellinger", distance="bray"){
  if (is.null(transformation)){
    dist <- vegdist(as.matrix(t(functionA(x))))
  } else {
    transform <- decostand(as.matrix(t(functionA(x))), transformation)
    dist <- vegdist(transform, method=distance)
  }
  pca <- prcomp(dist)
  pca
}

#' @export
betamatPCAvis <- function(x, geom="point", groups, reps, colours){
  pcaplot <- fviz_pca_ind(x, geom=geom, habeillage=rep(groups, reps),palette=colours, invisible="quali")+theme_bw()
  pcaplot
}

