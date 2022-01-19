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
fpath <- system.file("extdata", "toydata.xlsx", package="MetAlyzer")
```

```r
obj <- new("MetAlyzer")

obj <- init(obj, fpath, sheet = 1)
obj <- readData(obj)
obj <- readQuantStatus(obj)
obj <- filterMetabolites(obj)
show(obj)
summariseQuantData(obj)
## Not run: 
obj <- createPlottingData(obj, Method, Tissue)
obj <- imputePlottingData(obj, Tissue, Metabolite)
obj <- transformPlottingData(obj)
obj <- performANOVA(obj, "Method")
```



