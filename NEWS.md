# mirlyn 1.2.0

* Added a `mc.cores` option to `mirl()` (note that this change necessitated altering the behaviour of setting the random seed, so calls to the `mirl()` function will differ from 1.1.0 to 1.2.0)

* Fixed broken examples for `alphawichVis()`, `alphaconeVis()`

* General code, documentation cleanup

# mirlyn 1.1.0

* New function `get_asv_table()`: create an ASV table with abundance and taxonomy

* New function `phyloseq_to_df()`: create a `data.frame` containing abundance, taxonomy and metadata from `phyloseq` objects

* New function `plot_heat()`: create heatmaps visualizing relative abundance of a taxonmic group

* New function `randomseqsig()`: determine whether a taxonomic group of interest is over/under-represented
