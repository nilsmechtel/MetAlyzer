MetAlyzer
========

**An R Package to read and analyze MetIDQ&trade; output**

<!-- badges: start -->
[![R-CMD-check](https://github.com/andresenc/MetAlyzer/workflows/R-CMD-check/badge.svg)](https://github.com/andresenc/MetAlyzer/actions)
[![license](https://img.shields.io/badge/license-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)
<!-- badges: end -->

The package provides methods to read output files from the MetIDQ&trade; software into R. Metabolomics data is read and reformatted into an S4 object for convenient data handling, statistics and downstream analysis.

## Install

There is a version available on CRAN.

```r
install.packages("MetAlyzer")
```

## Quick start

![Overview](vignettes/MetAlyzer_overview.png)

The package takes metabolomic measurements and the quantification status (e.g. "Valid", "LOQ", "LOD") as ".xlsx" files generated from the MetIDQ&trade; software. Additionally, meta data for each sample can be provided for further analysis.

#### Set data path and read meta data
```r
fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
mpath <- system.file("extdata", "example_meta_data.rds", package = "MetAlyzer")
```

#### Create MetAlyzer object:
```r
obj <- MetAlyzerDataset(file_path = fpath)
```

#### Show MetAlyzer object:
```r
show(obj)
-------------------------------------
File name: example_data.xlsx 
Sheet: 1 
File path: /Library/Frameworks/R.framework/Versions/4.0/Resources/library/MetAlyzer/extdata 
Metabolites: 862 
Classes: 24 
Including metabolism indicators: TRUE 
Number of samples: 74 
Columns meta data: "Plate Bar Code"; "Sample Bar Code"; "Sample Type"; "Group"; "Tissue"; "Sample Volume"; "Measurement Time"
Ploting data created: FALSE
```

#### Use filter functions to exclude the metabolite indicators and only keep Group 1 to 6:
```r
obj <- filterMetabolites(obj, class_name = "Metabolism Indicators")
232 metabolites were filtered!
obj <- filterMetaData(obj, column = Group, keep = c(1:6))
```

#### Show statistics:
```r
summariseQuantData(obj)
-------------------------------------
Valid: 21951 (48.39%)
LOQ: 2832 (6.24%)
LOD: 20577 (45.36%)
NAs: 0 (0%)
```

#### Add meta data:
```r
meta_df <- readRDS(mpath)
obj <- updateMetaData(obj, "Replicate", meta_df$Replicate)
```

#### Reformat for plotting:
For further filtering and plotting, the data can be reformatted into a data frame.
```{r}
obj <- createPlottingData(obj, Group, Tissue)
gg_df <- plottingData(obj)

head(gg_df)
# A tibble: 6 × 12
# Groups:   Group, Tissue, Metabolite [2]
  Group Tissue     Metabolite Class          Replicates Concentration  Mean    SD    CV CV_thresh Status Valid
  <fct> <fct>      <fct>      <fct>               <int>         <dbl> <dbl> <dbl> <dbl> <fct>     <fct>  <lgl>
1 1     Drosophila C0         Acylcarnitines          1         203    179.  82.4 0.461 more30    Valid  TRUE 
2 1     Drosophila C0         Acylcarnitines          2          86.8  179.  82.4 0.461 more30    Valid  TRUE 
3 1     Drosophila C0         Acylcarnitines          3         246    179.  82.4 0.461 more30    Valid  TRUE 
4 2     Drosophila C0         Acylcarnitines          1         198    231. 124.  0.538 more30    Valid  TRUE 
5 2     Drosophila C0         Acylcarnitines          2         369    231. 124.  0.538 more30    Valid  TRUE 
6 2     Drosophila C0         Acylcarnitines          3         127    231. 124.  0.538 more30    Valid  TRUE 
```

#### Plot filter data and plot concentration of glutamic acid:
```{r}
glu_gg_df <- filter(gg_df, Metabolite == "Glu")

ggplot(glu_gg_df, aes(Group, Concentration, color = Status)) +
  geom_point() +
  scale_color_manual(values = c("Valid" = "#00CD66",
                                "LOQ" = "#87CEEB",
                                "LOD" = "#6A5ACD")) + 
  facet_grid(~ Tissue)
```
![](vignettes/example_ggplot.png)


## Detailed instructions
**For a comprehensive tutorial, please check out the vignette.**
