# Multiple Iterations of Rarefying for Library Normalization (mirlyn)

Differences in library sizes of amplicon sequencing datasets do not represent true biological variation. The differences observed between library sizes of different samples requires a method to control and correct for to allow for diversity analyses to be conducted accurately. Rarefying is a common normalization technique that implements the random subsampling of libraries to create rarefied libraries.
However, despite the frequent usage of rarefying, it has recently been criticized for being statistically inadmissable due to the omission of valid data (McMurdie and Holmes, 2014). While we acknowledge that rarefying has the potential to introduce variation through the omission of sequences, rarefying is a statistical tool that when handled carefully and with the necessary understanding of your data
that can be implemented succesfully for diversity analyses when applied over multiple iterations. This package, mirlyn, provides the necessary functions to conduct diversity analyses with repeated multiple iterations of rarefying to characterize the uncertainty introduced through the random subsampling of rarefying.

## Getting Started

A variety of normalization techniques have been proposed for amplicon sequencing data. It is strongly advised that you review literature resources to ensure that you make an educated decision about the normalization techniques you are applying to your data. 

### Suggested Resources
McMurdie and Holmes
Gloor et al
Weiss et al
Shade et al

## Functions

mirlyn includes a suite of functions focusing on the application of multiple iterations of rarefying for library normalization but also includes other functions that are useful for taxonomic marker gene community analysis. 
- asv_rename() - Assigns unique ASV I.D's to sequence variants.
- filt_fasta_rename() - Assigns corresponding unique ASV identifiers to a FASTA file.
- bartax() - Function for generating taxonomic composition bar charts for specified taxonomic levels. 
- fullbartax() - Function for generating taxonomic composition bar charts for all levels of taxonomy. 
- rarecurve() - Generates a rarefaction curve for sample sets. 
- mirl() - Repeatedly rarefies samples. 
- alphawichDF() - Generates a dataframe containing specified alpha diversity indices from mirl() output.
- alphawichVis() - Generates a ggplot2 object visualizing the alpha-diversity index from dataframe. 
- alphacone() - Generates a visualization of the distribution of diversity index for all different library sizes. 
- betamatPCA() - Generates a PCA for mirl() output.
- betamatPCAvis() - Generates a visualization for PCA output.
- divana() - All in one function to conduct taxonomic composition bar charts and full diversity analyses on samples. 



## Installation

mirlyn is available via github:

```
devtools:install_github("escamero/mirlyn")
```

This is a test sentence.
