MetAlyzer
========

**An R Package to read and analyze MetIDQ:tm: output**

The package provides functions to read output files from the MetIDQ:tm:software into R. Metabolomics and meta data are read and reformatted into an S4 object for convenient data handling, statistics and downstream analysis.

## Install

There is a version available on CRAN.

```r
install.packages('MetAlyzer')
```

## Quick start

![Overview](vignettes/MetAlyzer_overview.png)

The package takes ".xlsx" files generated as output from the MetIDQ:tm:software. Additionally, meta data for each sample can be provided for further analysis.

Set data path and read meta data:
```r
fpath <- system.file("data", "example.xlsx", package="MetAlyzer")
mpath <- system.file("data", "example_meta_data.rds", package="MetAlyzer")
```

Initialize MetAlyzer object and load data
```r
obj <- new("MetAlyzer")
obj <- init(obj, file_path=fpath)
obj <- readData(obj)
```

Show MetAlyzer object
```r
> show(obj)
-------------------------------------
File name: example.xlsx 
Sheet: 1 
# File path: ~/MetAlyzer/data 
Metabolites: 40 
Classes: 1 
Including metabolism indicators: FALSE 
Number of samples: 18 
Columns meta data: "Plate Bar Code"; "Sample Bar Code"; "Sample Type"; "Group"; "Sample Volume"; "Measurement Time"
Ploting data created: FALSE 
```

Show statistics
```r
> summariseQuantData(obj)
-------------------------------------
Valid: 248 (34.44%)
LOD: 423 (58.75%)
LOQ: 49 (6.81%)
NAs: 0 (0%)
```

Add meta data
```r
obj <- updateMetaData(obj, "Replicate", as.vector(meta_df$Replicate))
```



