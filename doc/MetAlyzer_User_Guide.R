## ---- echo=FALSE, message=FALSE-----------------------------------------------
knitr::opts_chunk$set(echo=TRUE, collapse=T, comment='#>')
library(MetAlyzer)
library(ggplot2)
library(dplyr)

## ----set_data_path------------------------------------------------------------
fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
mpath <- system.file("extdata", "example_meta_data.rds", package = "MetAlyzer")

## ----initialize---------------------------------------------------------------
obj <- MetAlyzerDataset(file_path = fpath)
show(obj)

## ----getter-------------------------------------------------------------------
head(metaData(obj))

head(quantStatus(obj), c(5, 5))

## ----filter_metabolites-------------------------------------------------------
obj <- filterMetabolites(obj, class_name = "Metabolism Indicators")

## ----filter_meta_data---------------------------------------------------------
obj <- filterMetaData(obj, column = Group, keep = c(1:6))

## ----summarise_quant_data-----------------------------------------------------
summariseQuantData(obj)

## ----read_meta_data-----------------------------------------------------------
meta_df <- readRDS(mpath)
head(meta_df)

## ----update_meta_data---------------------------------------------------------
obj <- updateMetaData(obj, name = Replicate, new_colum = meta_df$Replicate)
show(obj)

## ----rename_meta_data---------------------------------------------------------
obj <- renameMetaData(obj, Method = Group)
head(metaData(obj))

## ----plotting_data------------------------------------------------------------
obj <- createPlottingData(obj, Method, Tissue, ungrouped = Replicate)
gg_df <- plottingData(obj)

head(gg_df)

## ---- glu_plot, fig.width=7, fig.height=4.5-----------------------------------
glu_gg_df <- filter(gg_df, Metabolite == "Glu")

ggplot(glu_gg_df, aes(Method, Concentration, color = Status)) +
  geom_point() +
  scale_color_manual(values = c("Valid" = "#00CD66",
                                "LOQ" = "#87CEEB",
                                "LOD" = "#6A5ACD")) + 
  ylab("Concentration [pmol/mg Tissue]") + 
  facet_grid(~ Tissue)

## ----before_imp---------------------------------------------------------------
length(which(gg_df$Concentration == 0))

## ----imputation---------------------------------------------------------------
obj <- imputePlottingData(obj, Tissue, Metabolite)
imp_gg_df <- plottingData(obj)

## ----after_imp----------------------------------------------------------------
min(imp_gg_df$imp_Conc, na.rm = TRUE)

## ----transformation-----------------------------------------------------------
obj <- transformPlottingData(obj)
trans_gg_df <- plottingData(obj)
head(trans_gg_df)

## ---- glu_plot_transformed, fig.width=7, fig.height=4.5-----------------------
trans_glu_gg_df <- filter(trans_gg_df, Metabolite == "Glu")

ggplot(trans_glu_gg_df, aes(Method, transf_Conc, color = Status)) +
  geom_point() +
  scale_color_manual(values = c("Valid" = "#00CD66",
                                "LOQ" = "#87CEEB",
                                "LOD" = "#6A5ACD")) +
  facet_grid(~ Tissue)

## ---- ANOVA-------------------------------------------------------------------
obj <- performANOVA(obj, categorical = Method)
anv_gg_df <- plottingData(obj)
head(data.frame(anv_gg_df))

## ---- ANOVA_optimal-----------------------------------------------------------
anv_gg_df$optimal <- sapply(anv_gg_df$ANOVA_group, function(g) grepl("A", g))
obj <- setPlottingData(obj, anv_gg_df)
head(data.frame(anv_gg_df))

