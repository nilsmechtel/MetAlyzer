MetAlyzer
========

** An R Package to read and analyze MetIDQ:tm: output**

The package provides R functions to read output files from the MetIDQ:tm: software. 

## Install

There is a version available on CRAN.

```r
install.packages('MetAlyzer')
```

## Quick start

![Overview](`r rprojroot::is_r_package$find_file('inst/images/MetAlyzer_overview.png')`)

The 

```r
obj <- new("MetAlyzer")
fpath <- system.file("extdata", "toydata.xlsx", package="MetAlyzer")
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



