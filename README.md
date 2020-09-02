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

## Installation

mirlyn is available via github:

```
devtools:install_github("escamero/mirlyn")
```

## Functions

mirlyn includes a suite of functions focusing on the application of multiple iterations of rarefying for library normalization but also includes other functions that are useful for taxonomic marker gene community analysis. 
- asv_rename() - Assigns unique ASV I.D's to sequence variants.
- fasta_rename() - Assigns corresponding unique ASV identifiers to a FASTA file.
- bartax() - Function for generating taxonomic composition bar charts for specified taxonomic levels. 
- fullbartax() - Function for generating taxonomic composition bar charts for all levels of taxonomy. 
- rarecurve() - Generates a rarefaction curve for sample sets. 
- mirl() - Repeatedly rarefies samples. 
- alphawichDF() - Generates a dataframe containing specified alpha diversity indices from mirl() output.
- alphawichVis() - Generates a ggplot2 object visualizing the alpha-diversity index from dataframe. 
- alphacone() - Generates a visualization of the distribution of diversity index for all different library sizes. 
- betamatPCA() - Generates a PCA for mirl() output.
- betamatPCAvis() - Generates a visualization for PCA output.

## Example Workflow

### 1. Assigning ASV Identifiers to Sequence Variants

Prior to conducting analysis of sequencing data, users may choose to assign unique and easily referenced identifiers to their data. 

```
library(mirlyn)
library(phyloseq)
data(GlobalPatterns)
asv_rename(GlobalPatterns)
```

After unique and easily referenced ASV identifiers are assigned to sequencing data, users can cross-reference the original sequence variant identifiers in the FASTA file to the new identifiers using the fasta_rename function.


### 2. Taxonomic Composition Bar Charts

Prior to conducting diversity analysis, users may choose to generate relative abundance taxonomic composition bar charts. While this is not the focus of mirlyn, visualizations are supported. Users may choose to generate compositional bar charts for one taxonomic level or may generate graphs for all taxonomic levels. Since taxonomic bar charts utilize relative abundance, no rarefaction is required. Users may choose to specify the colours to use which must be equal to the number of taxonomic groups included in the generated graph. If users choose to generate all taxonomic graphs at once, graphs will be generated in a list object. 

```
#Generating taxonomic bar graph for Phylum taxonomic rank.
library(mirlyn)
data(example)
cols <- c("black", "darkgoldenrod1", "dodgerblue", "deeppink4", "chartreuse3", "burlywood4", "navy", "blueviolet", "tan2", "lavenderblush3", "cyan4")

bartax(example, "Sample", taxrank = "Phylum", cols = cols)
```

```
#Generating taxonomic bar graphs for all taxonomic ranks.
library(mirlyn)
data(example)
fullbartax(example, "sample")
Class
```

### 3. Generate Rarefaction Curve to Determine Appropriate Rarefied Library Size

Prior to rarefying, a rarefaction curve can be generated to provide users with an overview of the observed sequence variants in samples corresponding to different rarefied library sizes. Theoretically, samples that have a plateau in the rarefaction curve have maximal observed diversity. Rarefaction curves should be used in rarefied library size selection to ensure that the user-selected library size encompasses maximal diversity while being inclusive of all samples where possible. To generate the rarefaction curves for a sample set, users must first rarefy the entire sample set followed by generation of the rarefaction curve. 

```
#Rarefy entire sample set
library(mirlyn)
data(example)
rarefy_whole_rep_ex <- rarefy_whole_rep(example, rep = 100)

#Generate rarefaction curve
rarecurve_ex <- rarecurve(rarefy_whole_rep_ex, sample = "Sample")
```

### 4. Multiple Iterations of Rarefying Libraries (mirl)

After generating rarefaction curves, users may select an appropriate rarefied library size for their analysis. Users should aim to select a library size that represents maximal diversity and is inclusive of all samples. In the case where users must make the decision between losing samples or drastically reducing the represented diversity, users may opt to conduct analyses at the lower library size inclusive of all samples at the loss of diversity in some samples in addition to a larger rarefied library size which results in exclusion of small library size samples. Depending on the data structure, users may choose to include a different number of repeated iterations. For example, if the repeated iterations do not result in highly variable outputs in the diversity analyses, the number of iterations may be reduced. However, if large variation is present, users should aim to include a larger number of iterations to allow for better characterization of variation introduced through random subsampling. The mirl_object will be used in the subsequent analyses. 

```
library(mirlyn)
data(example)
mirl_object <- mirl(example, libsize = 10000, rep = 100, set.seed = 120)
```

### 5. Alpha-Diversity

mirlyn contains two visualization options for alpha-diversity analyses. Both implement the use of a diversity metric (e.g., Shannon diversity index). The alphawichDF() function utilizes the mirl_object generated in the previous step and is only applicable to the diversity metric at the specified library size used with mirl(). The alphacone() function generates a distrubtion of the diversity metric across different rarefied library sizes providing users with a comprehensive view of the diversity metric as a function of rarefied library size. 

```
#Alphawich 
alphadiv_df <- alphadivDF(mirl_object)
alphawichVis(example, "Sample")

#Alphacone
data(example)
alphacone_example <- alphacone(example, rep = 100)
alphaconeVis(alphacone_example, "Sample")
```

### 6. Beta-Diversity

Currently, mirlyn only supports the use of PCA for beta-diversity analyses. Future ordination techniques such as PCoA and NMDS may be implemented in future versions. A Hellinger transformation is recommended to apply to sequence count data prior to conducting PCA to account for the arch-effect regularly seen in ecological data. The beta-diversity functions utilize the mirl_object generated previously. 

```
betamatPCA_object <- betamatPCA(mirl_object, dsim = "bray")
betamatPCAvis(betamatPCA_object, groups = c("A", "B", "C", "D", "E", "F), reps = 100, colours = c("#000000", "#E69F00", "#0072B2", "#009E73", "#F0E442", "#D55E00"))
```

