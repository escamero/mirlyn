% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomseqsig.R
\name{randomseqsig}
\alias{randomseqsig}
\title{Compositional Significance Testing}
\usage{
randomseqsig(x, taxlevel = "Phylum", group = "Cyanobacteria", nshuff = 1000)
}
\arguments{
\item{x}{is a \code{phyloseq} object.}

\item{taxlevel}{is the taxonomic level of the group of interest.}

\item{group}{is the name of the taxonomic group of interest.}

\item{nshuff}{is the number of times data shuffling is performed.}
}
\value{
A \code{list} object which contatins: 1) \code{objshuff}: The shuffled sequence
variant counts. 2) \code{stats}: The mean, standard deviation, real counts and
p-values of the shuffled data.
}
\description{
The \code{\link[=randomseqsig]{randomseqsig()}} function will identify whether a taxonomic group of interest
is significantly over-represented or underrepresented in the community.
}
\examples{
library(mirlyn)
data(example)
significance_example <- randomseqsig(example, nshuff = 10)

}
