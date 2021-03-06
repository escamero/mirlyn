% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bartax.R
\name{bartax}
\alias{bartax}
\title{Taxonomic Composition Stacked Barcharts}
\usage{
bartax(
  x,
  xvar = "Sample",
  lowabun = 0.01,
  yvar = "Abundance",
  taxrank = "Genus",
  cols = NULL
)
}
\arguments{
\item{x}{is a \code{phyloseq} object.}

\item{xvar}{is a categorical variable selected from the metadata.}

\item{lowabun}{is the cut off value for a low-abundant group. Set to \code{NULL} to not use.}

\item{yvar}{is the relative abundance of reads.}

\item{taxrank}{is the specified taxonomic rank.}

\item{cols}{is a vector containing specified colour values (sufficient number for
number of taxonomic groups).}
}
\value{
A \code{ggplot} object.
}
\description{
The \code{\link[=bartax]{bartax()}} function will create taxonomic composition barcharts. If desired,
low-abundant taxonomic groups will be clustered into a single low abundant group.
}
\examples{
library(mirlyn)
data(example)
cols <- c("black","darkgoldenrod1","dodgerblue","deeppink4","chartreuse3",
   "burlywood4","navy","blueviolet", "tan2","lavenderblush3","cyan4")
bartax(example, "Sample", lowabun = 0.01, yvar = "Abundance", taxrank = "Phylum", cols)

}
