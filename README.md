# Update: Version 1.3.0 Now Available

# Multiple Iterations of Rarefying for Library Normalization (mirlyn)

Differences in library sizes of amplicon sequencing datasets do not represent true biological variation. The differences observed between library sizes of different samples requires a method to control and correct for to allow for diversity analyses to be conducted accurately. Rarefying is a common normalization technique that implements the random subsampling of libraries to create rarefied libraries.
However, despite the frequent usage of rarefying, it has recently been criticized for being statistically inadmissable due to the omission of valid data (McMurdie and Holmes, 2014). While we acknowledge that rarefying has the potential to introduce variation through the omission of sequences, rarefying is a statistical tool that when handled carefully and with the necessary understanding of your data
that can be implemented succesfully for diversity analyses when applied over multiple iterations. This package, mirlyn, provides the necessary functions to conduct diversity analyses with repeated multiple iterations of rarefying to characterize the uncertainty introduced through the random subsampling of rarefying.

## Getting Started

A variety of normalization techniques have been proposed for amplicon sequencing data. It is strongly advised that you review literature resources to ensure that you make an educated decision about the normalization techniques you are applying to your data. 

### Suggested Resources
Cameron, E.S., Schmidt, P.J., Tremblay, B.J.-M., Emelko, M.B., Müller, K.M., 2020. To rarefy or not to rarefy: Enhancing microbial community analysis through next-generation sequencing. bioRxiv, https://doi.org/10.1101/2020.09.09.290049


Gloor, G.B., Macklaim, J.M., Pawlowsky-Glahn, V., Egozcue, J.J., 2017. Microbiome datasets are compositional: And this is not optional. Front. Microbiol. 8, 1–6. https://doi.org/10.3389/fmicb.2017.02224


McMurdie, P.J., Holmes, S., 2014. Waste not, want not: Why rarefying microbiome data is inadmissible. PLoS Comput. Biol. 10. https://doi.org/10.1371/journal.pcbi.1003531


Weiss, S., Xu, Z.Z., Peddada, S., Amir, A., Bittinger, K., Gonzalez, A., Lozupone, C., Zaneveld, J.R., Vázquez-Baeza, Y., Birmingham, A., Hyde, E.R., Knight, R., 2017. Normalization and microbial differential abundance strategies depend upon data characteristics. Microbiome 5, 1–18. https://doi.org/10.1186/s40168-017-0237-y

## Installation

mirlyn is available via github:

```r
install.packages("BiocManager")
BiocManager::install("escamero/mirlyn")
```

## Functions

mirlyn includes a suite of functions focusing on the application of multiple iterations of rarefying for library normalization but also includes other functions that are useful for taxonomic marker gene community analysis. 
- `asv_rename()` - Assigns unique ASV I.D's to sequence variants.
- `fasta_rename()` - Assigns corresponding unique ASV identifiers to a FASTA file.
- `bartax()` - Function for generating taxonomic composition bar charts for specified taxonomic levels. 
- `fullbartax()` - Function for generating taxonomic composition bar charts for all levels of taxonomy. 
- `rarecurve()` - Generates a rarefaction curve for sample sets. 
- `mirl()` - Repeatedly rarefies samples. 
- `alphadivDF()` - Generates a dataframe containing specified alpha diversity indices from `mirl()` output.
- `alphawichVis()` - Generates a ggplot2 object visualizing the alpha-diversity index from dataframe. 
- `alphacone()` - Generates a visualization of the distribution of diversity index for all different library sizes. 
- `betamatPCA()` - Generates a PCA for `mirl()` output.
- `betamatPCAvis()` - Generates a visualization for PCA output.
- `phyloseqtodf()` - Generations a dataframe from a `phyloseq` object including taxonomy, ASV read counts and supporting metadata. 
- `get_asv_table()` - Generates a compiled ASV table with counts and taxonomic classification of individual amplicon sequence variants. 
- `randomseqsig()` - Identification of whether a taxonomic group of interest is significantly dominant or rare in the community using data shuffling. 
- `plot_heat()` - Generates heat maps visualizing the relative abundance of taxonomic group of interest from a dataframe.

To view help page for any of the functions type the function name preceded by a question mark. For example: `?mirl`

## Example R Workflow

## a) Data Import
`QIIME2R` was used to import .qza files into R as a `phylsoeq` objects. This generated phyloseq object was then used for community diversity analyses using `mirlyn`, an R package developed in Chapter 2 including various functions for library normalization and diversity analyses. `mirlyn` funcitons frequently relies on the usage of a `phyloseq` object and various methods can be used for generating this object in R. 

## b) Data Handling
### 1. Assigning ASV Identifiers to Sequence Variants

During the creation of the ASV table in QIIME2, each unique sequence is assigned an identifier consisting of a string of characters (e.g., 88b44c11059bcf2950ca0ac50f3eb08f). To improve readability, the asv_rename() function codes these character string identifiers to a new identifier in the form “ASV###”. While this step is not mandatory, it allows for easy reference to specific ASV. 
```r
library(mirlyn)
library(phyloseq)

# Load example dataset from phyloseq.
data(GlobalPatterns) 

# Renames ASV in phyloseq object
asv_rename(GlobalPatterns)
```

After unique and easily referenced ASV identifiers are assigned to sequencing data, users can cross-reference the original sequence variant identifiers in the FASTA file to the new identifiers using the fasta_rename function.

### 2. Generation of Data Frame from `phyloseq` Object
Data is initially imported into R as a `phyloseq` object. The `phyloseq` object is critical for subsequently performing diversity analysis. However, for plotting options or subsequent export as a CSV file, the `phyloseq_to_df()` function will convert the phyloseq object to a data frame containing the ASV counts, taxonomic classification and metadata.  

```r
# Generate dataframe from phyloseq object
example_df <- phyloseq_to_df(example)
```

### 3. Generation of Compiled ASV Table
The data frame generated using `phyloseq_to_df()` can be further organized to focus on the read counts of each ASV across the different samples. The `get_asv_table()` also includes the taxonomic classification of the ASV but does not include sample metadata. 

```r
# Generate compiled ASV table from phyloseq dataframe
example_df_asv <- get_asv_table(example_df)
```

## c) Taxonomic Composition
### 1. Visualization

`mirlyn` provides two visualization options for taxonomic communities including stacked bar plots and heatmaps. Heatmaps are optimal to use when interested in exploring the trends in relative abundances of one taxonomic group (e.g., Cyanobacteria). Alternatively, the stacked bar charts can be used to identify overall composition of communities. 

```r
# Stacked Barcharts at the Phylum Rank
cols <- c("black", "darkgoldenrod1", "dodgerblue", "deeppink4", "chartreuse3", "burlywood4", "navy", "blueviolet", "tan2", "lavenderblush3", "cyan4")

example_barchart <- bartax(example, “Sample”, taxrank = “Phylum”, cols = cols) 

# Heatmap of Cyanobacterial Abundances in the Bacterial Community
example_df_phylum <- example_df %>% group_by(Sample, Id, Phylum) %>% summarise(abaundance = (sum(abundance)) %>% mutate(Proportion = abundance/sum(abundance)*100)

plot_heat(example_df_phylum, taxlevel = “Phylum”, taxaname = “Cyanobacteria”, xvar = “sample”, yvar = “Id”, fillvar = “Proportion”)+scae_fill_gradient(low = “white”, high = “midnightblue”)
```

```r
# Generating taxonomic bar graphs for all taxonomic ranks.
library(mirlyn)

# Load example data from mirlyn.
data(example)

# Generates taxonomic barchart for all taxonomic levels as a list object. 
alltaxgraphs <- fullbartax(example, "Sample")

# To show all availabe graphs
names(alltaxgraphs)

# To look at the graph for a specific taxonomic level, for example: Class. 
alltaxgraphs$Class
```

### 2. Compositional Significance Testing
Amplicon sequencing data is inherently compositional (Gloor et al., 2016). The composition of these communities is reported in relative abundance but raises the question of when is a group statistically abundant within the community. The `randomseqsig()` function will identify whether a taxonomic group of interest is significantly dominant in the community. This can be used to identify conditions where a taxonomic group of interest (e.g., Cyanobacteria) are in significantly higher abundances. 
```r
# Calculate significance of Phylum: Cyanobacteria
compsig_example <- randomseqsig(example, taxlevel = "Phylum", group = "Cyanobacteria", nshuff = 1000)
```

## d) Diversity Analyses
Prior to conducting diversity analysis, libraries must be normalized to account for variation in library sizes. A variety of techniques are available each with their own benefits and limitations and researchers are encouraged to evaluate the effectiveness of these techniques for their data. However, for this research, `mirlyn` utilizes repeated iterations of rarefying, the process of subsampling to a user specified library size.

### 1. Library Normalization
To identify appropriate library sizes to rarefy to, a rarefaction curve can be generated to provide an overview of the observed ASV in samples corresponding to different rarefied library sizes. Theoretically, samples that have a plateau in the curve have reached maximal observed diversity. This visualization should be used to select an appropriate library size which encompasses maximal diversity while being inclusive of samples. 

```r
# Creation of rarefaction curve data frame
Rarefy_whole_rep_example <- rarefy_whole_rep(example, rep = 100)

#Visualization of rarefaction curve
Rarecurve_ex <- rarecurve(rarefy_whole_rep_ex, sample = “Sample”
```

### 2. Multiple Iterations of Rarefying Libraries (mirl)

After generating rarefaction curves, users may select an appropriate rarefied library size for their analysis. Users should aim to select a library size that represents maximal diversity and is inclusive of all samples. In the case where users must make the decision between losing samples or drastically reducing the represented diversity, users may opt to conduct analyses at the lower library size inclusive of all samples at the loss of diversity in some samples in addition to a larger rarefied library size which results in exclusion of small library size samples. Depending on the data structure, users may choose to include a different number of repeated iterations. For example, if the repeated iterations do not result in highly variable outputs in the diversity analyses, the number of iterations may be reduced. However, if large variation is present, users should aim to include a larger number of iterations to allow for better characterization of variation introduced through random subsampling. The `mirl_object` will be used in the subsequent analyses. 

```r
library(mirlyn)

# Load example data from mirlyn.
data(example)

# Repeatedly rarefies to specified library size n times
mirl_object <- mirl(example, libsize = 10000, rep = 100, set.seed = 120)
```

### 3. Alpha-Diversity

mirlyn contains two visualization options for alpha-diversity analyses. Both implement the use of a diversity metric (e.g., Shannon diversity index). The `alphadivDF()` function utilizes the `mirl_object` generated in the previous step and is only applicable to the diversity metric at the specified library size used with `mirl()`. The `alphacone()` function generates a distrubtion of the diversity metric across different rarefied library sizes providing users with a comprehensive view of the diversity metric as a function of rarefied library size. 

```r
# Alphawich Functions
# Generates dataframe of alpha-diversity metric from mirl_object
alphadiv_df <- alphadivDF(mirl_object)

# Generates visualization from alphadiv_df. Substitute xvar for your own metadata column. 
alphawichVis(alphadiv_df, xvar = "Id")

#Alphacone Functions
# Load example data from mirlyn.
data(example)

# Generates dataframe of alpha-diversity metric across all library sizes.
alphacone_example <- alphacone(example, rep = 100)

# Generates distribution plot of alpha-diversity metric across library sizes. 
alphaconeVis(alphacone_example, "Sample")
```

### 4. Beta-Diversity

Currently, mirlyn only supports the use of PCA for beta-diversity analyses. Future ordination techniques such as PCoA and NMDS may be implemented in future versions. A Hellinger transformation is recommended to apply to sequence count data prior to conducting PCA to account for the arch-effect regularly seen in ecological data. The beta-diversity functions utilize the `mirl_object` generated previously. 

```r
# Generates PCA from mirl_object for specified dissimilarity metric (dsim).
betamatPCA_object <- betamatPCA(mirl_object, dsim = "bray")

# Generates visualization of PCA.
betamatPCAvis(betamatPCA_object, groups = c("A", "B", "C", "D", "E", "F"), reps = 100, colours = c("#000000", "#E69F00", "#0072B2", "#009E73", "#F0E442", "#D55E00"))
```

