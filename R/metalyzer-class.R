library(xlsx)
library(dplyr)
library(tidyr)

#' A S4 class to read and analyze MetIDQ output
#'
#' @slot file_path A length-one character vector
#' @slot sheet A length-one numeric vector
#' @slot metabolites A character vector with all 630 measured metabolites and
#' optional 234 additional metabolism indicators
#' @slot raw_data A data frame containing all raw measurements;
#' dimension: # samples x # metabolites
#' @slot quant_status A data frame containing the quantification status of all
#' measurements; dimension: # samples x # metabolites
#' @slot meta_data A data frame containing any meta data;
#' dimension: # samples x # meta variables
#' @slot plotting_data A tibble data.frame containing reshaped information of raw_data,
#' quant_status and meta_data for plotting with ggplot2
#' @slot .full_sheet A matrix containing the un-sliced Excel sheet
#' @slot .data_ranges A length-six numeric list with rows and columns
#' information for slicing
#'
#' @return
#' @export
#'
#' @examples
#' obj <- new("MetAlyzer")
#' obj <- init(obj, file_path="results.xlsx")
#' obj <- readData(obj)
#' obj <- readQuantStatus(obj)
#' obj <- filterMetabolites(obj)
#' show(obj)
#' summariseQuantData(obj)
#' obj <- createPlottingData(obj, Method, Tissue)
#' obj_organisms <- imputePlottingData(obj_organisms, Tissue, Metabolite)
#' obj_organisms <- transformPlottingData(obj_organisms)
#' obj_organisms <- performANOVA(obj_organisms, "Method")

setClass("MetAlyzer",
         slots=list(
           ## input ##
           file_path="character",
           sheet="numeric",
           ## output ##
           metabolites="character",
           raw_data="data.frame",
           quant_status="data.frame",
           meta_data="data.frame",
           plotting_data="data.frame",
           ## internal ##
           .full_sheet="matrix",
           .data_ranges="list",
           .orig_metabolites="character"
         )
)


#' Initialize
#'
#' This method initializes file path and sheet number
#' @param MetAlyzer MetAlyzer object
#' @param file_path A length-one character vector with the file path
#' @param sheet A length-one numeric vector with the sheet number
#'
#' @return
#' @export

setGeneric("init",
           function(MetAlyzer, file_path, sheet=1)
             standardGeneric("init")
           )
setMethod("init",
          "MetAlyzer",
          function(MetAlyzer, file_path, sheet) {
            MetAlyzer@file_path <- file_path
            MetAlyzer@sheet <- sheet
            return(MetAlyzer)
          })


#' Show a MetAlyzer object
#'
#' This method calls show_obj() and shows the MetAlyzer object
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

# setGeneric("show",
#            function(MetAlyzer)
#              standardGeneric("show")
# )
setMethod("show",
          "MetAlyzer",
          function(MetAlyzer) {
            show_obj(MetAlyzer)
          }
)


#' Open file and read data
#'
#' This method opens the given MetIDQ output Excel file and extracts metabolites,
#' raw data and meta data
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("readData",
           function(MetAlyzer)
             standardGeneric("readData")
           )
setMethod("readData",
          "MetAlyzer",
          function(MetAlyzer) {
            MetAlyzer@.full_sheet <- open_file(MetAlyzer)
            MetAlyzer@.data_ranges <- get_data_range(MetAlyzer)
            MetAlyzer@.orig_metabolites <- read_metabolties(MetAlyzer)
            MetAlyzer@metabolites <- MetAlyzer@.orig_metabolites
            MetAlyzer@raw_data <- read_raw_data(MetAlyzer)
            MetAlyzer@meta_data <- read_meta_data(MetAlyzer)
            return(MetAlyzer)
          }
)


#' Read quantification status
#'
#' This method gets the background color of each cell of the opened Excel file
#' and assigns it to the corresponding quantification status
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("readQuantStatus",
           function(MetAlyzer)
             standardGeneric("readQuantStatus")
           )
setMethod("readQuantStatus",
          "MetAlyzer",
          function(MetAlyzer) {
            df_BG <- read_BG_color(MetAlyzer)
            MetAlyzer@quant_status <- assign_quant_status(df_BG)
            return(MetAlyzer)
          }
)


#' Filter metabolites
#'
#' This method calls filter_metabolites() and filters metabolites
#' @param MetAlyzer MetAlyzer object
#' @param class_name A class to be filtered out
#'
#' @return
#' @export

setGeneric("filterMetabolites",
           function(MetAlyzer, class_name="Metabolism Indicators")
             standardGeneric("filterMetabolites")
)
setMethod("filterMetabolites",
          "MetAlyzer",
          function(MetAlyzer, class_name) {
            filter_metabolites(MetAlyzer, class_name)
          }
)


#' Reset metabolites
#'
#' This method resets the filtered metabolites
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("resetMetabolites",
           function(MetAlyzer)
             standardGeneric("resetMetabolites")
)
setMethod("resetMetabolites",
          "MetAlyzer",
          function(MetAlyzer) {
            MetAlyzer@metabolites <- MetAlyzer@.orig_metabolites
            return(MetAlyzer)
          }
)


#' Filter meta data
#'
#' This method calls filter_meta_data() and filters meta_data
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("filterMetaData",
           function(MetAlyzer, column, keep=NULL, remove=NULL)
             standardGeneric("filterMetaData")
)
setMethod("filterMetaData",
          "MetAlyzer",
          function(MetAlyzer, column, keep, remove) {
            filter_meta_data(MetAlyzer, column, keep, remove)
          }
)


#' Reset meta data
#'
#' This method resets the filter of meta_data
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("resetMetaData",
           function(MetAlyzer)
             standardGeneric("resetMetaData")
)
setMethod("resetMetaData",
          "MetAlyzer",
          function(MetAlyzer) {
            MetAlyzer@meta_data$Filter <- TRUE
            return(MetAlyzer)
          }
)


#' Update meta data
#'
#' This method adds another column to unfiltered meta_data
#' @param MetAlyzer MetAlyzer object
#' @param name A length-one character vector giving the new column name
#' @param new_colum A vector for the new column (length has to be same as the
#' number of samples)
#'
#' @return
#' @export

setGeneric("updateMetaData",
           function(MetAlyzer, name, new_colum)
             standardGeneric("updateMetaData")
)
setMethod("updateMetaData",
          "MetAlyzer",
          function(MetAlyzer, name, new_colum) {
            MetAlyzer@meta_data[,name] <- factor(NA, levels = levels(new_colum))
            MetAlyzer@meta_data[,name][MetAlyzer@meta_data$Filter == TRUE] <- new_colum
            return(MetAlyzer)
          }
)


#' Rename meta data
#'
#' This method renames a column of meta_data using rename {dplyr}
#' @param MetAlyzer MetAlyzer object
#' @param ... Use new_name = old_name to rename selected variables
#'
#' @return
#' @export

setGeneric("renameMetaData",
           function(MetAlyzer, ...)
             standardGeneric("renameMetaData")
)
setMethod("renameMetaData",
          "MetAlyzer",
          function(MetAlyzer, ...) {
            MetAlyzer@meta_data <- rename(MetAlyzer@meta_data, ...)
            return(MetAlyzer)
          }
)


#' Get meta data
#'
#' This method returns the meta_data data frame with filtered samples (rows) and
#' all columns
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("metaData",
           function(MetAlyzer)
             standardGeneric("metaData")
)
setMethod("metaData",
          "MetAlyzer",
          function(MetAlyzer) {
            get_filtered_data(MetAlyzer, "meta")
          }
)


#' Get raw data
#'
#' This method returns the raw_data data frame with filtered samples (rows) and
#' metabolites (columns)
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("rawData",
           function(MetAlyzer)
             standardGeneric("rawData")
           )
setMethod("rawData",
          "MetAlyzer",
          function(MetAlyzer) {
            get_filtered_data(MetAlyzer, "data")
          }
)


#' Get quantification status
#'
#' This method returns the quant_status data frame with filtered samples (rows)
#' and metabolites (columns)
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("quantStatus",
           function(MetAlyzer)
             standardGeneric("quantStatus")
           )
setMethod("quantStatus",
          "MetAlyzer",
          function(MetAlyzer) {
            get_filtered_data(MetAlyzer, "quant")
          }
)


#' Get metabolites
#'
#' This method returns the filtered metabolites vector
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("metabolites",
           function(MetAlyzer)
             standardGeneric("metabolites")
)
setMethod("metabolites",
          "MetAlyzer",
          function(MetAlyzer) {
            MetAlyzer@metabolites
          }
)


#' Summarize quantification status
#'
#' This method calls sum_quant_data() and summarizes quantification status
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("summariseQuantData",
           function(MetAlyzer)
             standardGeneric("summariseQuantData")
           )
setMethod("summariseQuantData",
          "MetAlyzer",
          function(MetAlyzer) {
            sum_quant_data(MetAlyzer)
          }
)


#' Create plotting data
#'
#' This method combines raw_data, quant_status and meta_data to one tibble
#' data frame
#' @param MetAlyzer MetAlyzer object
#' @param ... A selection of columns from meta_data to add to reshaped data frame
#' @param ts A numeric vector of thresholds for CV categorization
#' @param valid_vec A character vector containing each quantification status that
#' is considered to be a valid measurement
#' @param t A numeric threshold to determine valid measurements
#'
#' @return
#' @export

setGeneric("createPlottingData",
           function(MetAlyzer, ...,
                    ts = c(0.1, 0.2, 0.3),
                    valid_vec = c("Valid", "LOQ"), t = 0.5)
             standardGeneric("createPlottingData")
)
setMethod("createPlottingData",
          "MetAlyzer",
          function(MetAlyzer, ..., ts, valid_vec, t) {
            plotting_data <- plotting_data(MetAlyzer, ...)
            plotting_data <- calc_CV(plotting_data, ts = ts)
            plotting_data <- valid_measurement(plotting_data, valid_vec, t)
            MetAlyzer@plotting_data <- plotting_data
            cat("-------------------------------------\n")
            cat("Returning plotting data\n")
            return(MetAlyzer)
          }
)


#' Get plotting data
#'
#' This method returns the plotting_data tibble data.frame
#' @param MetAlyzer MetAlyzer object
#'
#' @return
#' @export

setGeneric("plottingData",
           function(MetAlyzer)
             standardGeneric("plottingData")
)
setMethod("plottingData",
          "MetAlyzer",
          function(MetAlyzer) {
            MetAlyzer@plotting_data
          }
)


#' Impute plotting data
#'
#' This method adds the column imp_Conc to plotting_data containing imputed
#' concentration values (Concentration). Imputation: minimal positive value * i
#' @param MetAlyzer MetAlyzer object
#' @param i A numeric value below 1)
#'
#' @return
#' @export

setGeneric("imputePlottingData",
           function(MetAlyzer, ..., i = 0.2)
             standardGeneric("imputePlottingData")
)
setMethod("imputePlottingData",
          "MetAlyzer",
          function(MetAlyzer, ..., i) {
            plotting_data <- MetAlyzer@plotting_data
            grouping_vars <- group_vars(plotting_data)
            plotting_data <- plotting_data %>%
              group_by(...) %>%
              mutate(imp_Conc = zero_imputation(Concentration, i),
                     .after = Concentration) %>%
              group_by_at(grouping_vars)
            MetAlyzer@plotting_data <- plotting_data
            return(MetAlyzer)
          }
)


#' log transform plotting data
#'
#' This method adds the column transf_Conc to plotting_data with the logarithmic
#' transformation of imputed concentration values (imp_Conc)
#' @param MetAlyzer MetAlyzer object
#' @param func A logarithmic function to transform concentration values with
#' number of samples)
#'
#' @return
#' @export

setGeneric("transformPlottingData",
           function(MetAlyzer, func = log2)
             standardGeneric("transformPlottingData")
)
setMethod("transformPlottingData",
          "MetAlyzer",
          function(MetAlyzer, func) {
            MetAlyzer@plotting_data <- mutate(MetAlyzer@plotting_data,
                                              transf_Conc = log_transform(imp_Conc, func),
                                              .after = imp_Conc)
            return(MetAlyzer)
          }
)


#' ANOVA
#'
#' This method performs a one-way ANOVA adds the column Group to plotting_data
#' with the results of a Tukey post-hoc test
#' @param MetAlyzer MetAlyzer object
#' @param categorical A length-one character vector defining the categorical
#' variable)
#'
#' @return
#' @export

setGeneric("performANOVA",
           function(MetAlyzer, categorical)
             standardGeneric("performANOVA")
)
setMethod("performANOVA",
          "MetAlyzer",
          function(MetAlyzer, categorical) {
            plotting_data <- MetAlyzer@plotting_data
            grouping_vars <- group_vars(plotting_data)
            group_cols <- grouping_vars[-which(grouping_vars == categorical)]
            plotting_data <- plotting_data %>%
              group_by_at(group_cols) %>%
              mutate(Group = calc_anova(!!sym(categorical), transf_Conc, Valid)) %>%
              group_by_at(grouping_vars)
            MetAlyzer@plotting_data <- plotting_data
            return(MetAlyzer)
          }
)


#' Update plotting data
#'
#' This method adds or modifies a column to/of plotting_data
#' @param MetAlyzer MetAlyzer object
#' @param name A length-one character vector giving the new column name
#' @param new_colum A vector for the new column (length has to be same as the
#' number of samples)
#'
#' @return
#' @export

setGeneric("updatePlottingData",
           function(MetAlyzer, name, new_colum)
             standardGeneric("updatePlottingData")
)
setMethod("updatePlottingData",
          "MetAlyzer",
          function(MetAlyzer, name, new_colum) {
            MetAlyzer@plotting_data[,name] <- new_colum
            return(MetAlyzer)
          }
)
