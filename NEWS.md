# mirlyn 1.4.2

* Fix bug in `mirl()` introduced in v1.4.1

# mirlyn 1.4.1

* Custom prefix name in asvrenamer.R

* Don't forget user distance/similarity metric in `betamatPCA()`

* Warn about dropped samples in `mirl()`

# mirlyn 1.4.0

* Significant speed increase for internal function `repotu_df()`, which should help functions which use rarefied objects with many samples and reps

# mirlyn 1.3.1

* Fixed `plot_heat()` function

# mirlyn 1.3.0

* Speed ups to `alphadivDF()` and `betamatPCA()`

* Added set.seed and mc.cores options to `alphacone()`

# mirlyn 1.2.2

* Fixed another typo in README

* Added an informative error message to betamatPCAvis() when the number of colours and/or reps is incorrect

# mirlyn 1.2.1

* Fixed typos in README

* Fixed rng seed handling in `rarefy_whole_rep()` so that it is different when using multiple cores (note that this means the function is not reproducible between this version and previous ones)

# mirlyn 1.2.0

* Added a `mc.cores` option to `mirl()` (note that this change necessitated altering the behaviour of setting the random seed, so calls to the `mirl()` function will differ from 1.1.0 to 1.2.0)

* Fixed broken examples for `alphawichVis()`, `alphaconeVis()`

* General code, documentation cleanup

# mirlyn 1.1.0

* New function `get_asv_table()`: create an ASV table with abundance and taxonomy

* New function `phyloseq_to_df()`: create a `data.frame` containing abundance, taxonomy and metadata from `phyloseq` objects

* New function `plot_heat()`: create heatmaps visualizing relative abundance of a taxonmic group

* New function `randomseqsig()`: determine whether a taxonomic group of interest is over/under-represented
