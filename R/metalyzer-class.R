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
#' @slot plotting_data A tibble data frame containing reshaped information of raw_data,
#' quant_status and meta_data for plotting with ggplot2
#' @slot .full_sheet A matrix containing the un-sliced Excel sheet
#' @slot .data_ranges A length-six numeric list with rows and columns
#' information for slicing
#'
#' @export
#'
#' @examples
#' fpath <- system.file("extdata", "Extraction_test.xlsx", package="MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#' obj <- filterMetabolites(obj)
#' show(obj)
#' summariseQuantData(obj)
#'
#' obj <- renameMetaData(obj, Method = Group)
#' obj <- filterMetaData(obj, Method, keep = 1:6)
#'
#' obj <- createPlottingData(obj, Method, Tissue)
#' obj <- imputePlottingData(obj, Method, Metabolite)
#' obj <- transformPlottingData(obj)
#' \donttest{
#' obj <- performANOVA(obj, Methods)
#' }

MetAlyzer <- setClass("MetAlyzer",
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


#' Open file and read data
#'
#' This function creates a MetAlyzer object, opens the given MetIDQ output Excel
#' file and extracts metabolites, raw data, quantification status and meta data
#' @param file_path file path
#' @param sheet sheet index
#'
#' @return An MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- MetAlyzerDataset(file_path = "toydata.xlsx")
#' }

MetAlyzerDataset <- function(file_path, sheet=1) {
  MetAlyzer <- new("MetAlyzer", file_path = file_path, sheet = sheet)
  MetAlyzer@.full_sheet <- open_file(MetAlyzer)
  MetAlyzer@.data_ranges <- get_data_range(MetAlyzer)
  MetAlyzer@.orig_metabolites <- read_metabolties(MetAlyzer)
  MetAlyzer@metabolites <- MetAlyzer@.orig_metabolites
  MetAlyzer@raw_data <- read_raw_data(MetAlyzer)
  MetAlyzer@quant_status <- read_quant_status(MetAlyzer)
  MetAlyzer@meta_data <- read_meta_data(MetAlyzer)
  return(MetAlyzer)
}


#' Show a MetAlyzer object
#'
#' This method calls show_obj() and shows the MetAlyzer object
#' @param object MetAlyzer object
#'
#' @return A summary of the MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' show(obj)
#' # or
#' obj
#' }

setMethod("show",
          "MetAlyzer",
          function(object) {
            show_obj(object)
          }
)


#' Summarize quantification status
#'
#' This method calls sum_quant_data() and summarizes quantification status
#' @param MetAlyzer MetAlyzer object
#'
#' @return A summary of the quantification status
#'
#' @export
#'
#' @examples
#' \dontrun{
#' summariseQuantData(obj)
#' }

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


#' Filter metabolites
#'
#' This method calls filter_metabolites() and filters metabolites.
#' Note: metabo_vec overwrites class_name argument!
#' @param MetAlyzer MetAlyzer object
#' @param class_name A character value defining the class to be removed
#' @param metabo_vec A character vector defining metabolites to be removed
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- filterMetabolites(obj, class_name = "Metabolism Indicators")
#' obj <- filterMetabolites(obj, metabo_vec = c("C0", "C2", "C3"))
#' }

setGeneric("filterMetabolites",
           function(MetAlyzer, class_name="Metabolism Indicators",
                    metabo_vec=NULL)
             standardGeneric("filterMetabolites")
)
setMethod("filterMetabolites",
          "MetAlyzer",
          function(MetAlyzer, class_name, metabo_vec) {
            filter_metabolites(MetAlyzer, class_name, metabo_vec)
          }
)


#' Reset metabolites
#'
#' This method resets the filtered metabolites
#' @param MetAlyzer MetAlyzer object
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- resetMetabolites(obj)
#' }

setGeneric("resetMetabolites",
           function(MetAlyzer)
             standardGeneric("resetMetabolites")
)
setMethod("resetMetabolites",
          "MetAlyzer",
          function(MetAlyzer) {
            orig <- MetAlyzer@.orig_metabolites
            diff <- length(orig) - length(MetAlyzer@metabolites)
            if (diff > 0) {
              cat(paste("Restoring", diff, "metabolite(s).\n"))
            }
            MetAlyzer@metabolites <- orig
            return(MetAlyzer)
          }
)


#' Filter meta data
#'
#' This method calls filter_meta_data() and filters meta_data. Note: If both
#' "keep" and "remove" arguments are used "keep" overwrites the "remove"
#' argument.
#'
#' @param MetAlyzer MetAlyzer object
#' @param column A column of meta data for filtering
#' @param keep A vector defining which entries to keep from meta data
#' @param remove A vector defining which entries to remove meta data
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- filterMetaData(obj, Methods, keep = 1:6)
#' # or
#' obj <- filterMetaData(obj, Methods, remove = 7)
#' }

setGeneric("filterMetaData",
           function(MetAlyzer, column, keep=NULL, remove=NULL)
             standardGeneric("filterMetaData")
)
setMethod("filterMetaData",
          "MetAlyzer",
          function(MetAlyzer, column, keep, remove) {
            filter_meta_data(MetAlyzer, deparse(substitute(column)), keep, remove)
          }
)


#' Reset meta data
#'
#' This method resets the filter of meta_data
#' @param MetAlyzer MetAlyzer object
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- resetMetaData(obj)
#' }

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
#' @param name The new column name
#' @param new_colum A vector for the new column (length has to be same as the
#' number of filtered samples)
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- updateMetaData(obj, Date, format(Sys.Date()))
#' obj <- updateMetaData(obj, Analyzed, TRUE)
#' }

setGeneric("updateMetaData",
           function(MetAlyzer, name, new_colum)
             standardGeneric("updateMetaData")
)
setMethod("updateMetaData",
          "MetAlyzer",
          function(MetAlyzer, name, new_colum) {
            if (class(new_colum) == "factor") {
              levels <- levels(new_colum)
            } else {
              levels <- unique(new_colum)
            }
            meta_data <- MetAlyzer@meta_data
            chr_name <- deparse(substitute(name))
            meta_data[,chr_name] <- factor(NA, levels = levels)
            meta_data[,chr_name][meta_data$Filter == TRUE] <- new_colum
            MetAlyzer@meta_data <- meta_data
            return(MetAlyzer)
          }
)


#' Rename meta data
#'
#' This method renames a column of meta_data using rename {dplyr}
#' @param MetAlyzer MetAlyzer object
#' @param ... Use new_name = old_name to rename selected variables
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- renameMetaData(obj, Methods = X1)
#' }

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
#' @return The meta data data frame
#'
#' @export
#'
#' @examples
#' \dontrun{
#' metaData(obj)
#' }

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
#' @return The raw data data frame
#'
#' @export
#'
#' @examples
#' \dontrun{
#' rawData(obj)
#' }

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
#' @return The quantification status data frame
#'
#' @export
#'
#' @examples
#' \dontrun{
#' quantStatus(obj)
#' }

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
#' @return The metabolites vector
#'
#' @export
#'
#' @examples
#' \dontrun{
#' metabolites(obj)
#' }

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
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- createPlottingData(obj, Tissue, Method,
#' ts = c(0.1, 0.2, 0.3),
#' valid_vec = c("Valid", "LOQ"), t = 0.5)
#' }

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
            return(MetAlyzer)
          }
)


#' Get plotting data
#'
#' This method returns the plotting_data tibble data frame
#' @param MetAlyzer MetAlyzer object
#'
#' @return The plotting data data frame
#'
#' @export
#'
#' @examples
#' \dontrun{
#' plottingData(obj)
#' }

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
#' @param ... Variables to group by
#' @param i A numeric value below 1
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' imputePlottingData(obj, Tissue, Metabolite, i = 0.2)
#' }

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
#' number of samples
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- transformPlottingData(obj, func = log2)
#' obj <- transformPlottingData(obj, func = log)
#' }

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
#' @param categorical A  column defining the categorical variable
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- performANOVA(obj, Methods)
#' }

setGeneric("performANOVA",
           function(MetAlyzer, categorical)
             standardGeneric("performANOVA")
)
setMethod("performANOVA",
          "MetAlyzer",
          function(MetAlyzer, categorical) {
            plotting_data <- MetAlyzer@plotting_data
            grouping_vars <- group_vars(plotting_data)
            plotting_data <- plotting_data %>%
              group_by_at(grouping_vars, add = FALSE) %>%
              ungroup(Method) %>%
              mutate(ANOVA_group = calc_anova(Method, transf_Conc, Valid)) %>%
              group_by_at(grouping_vars)
            MetAlyzer@plotting_data <- plotting_data
            return(MetAlyzer)
          }
)


#' Update plotting data
#'
#' This method replaces plotting_data with an updated version
#' @param MetAlyzer MetAlyzer object
#' @param plotting_data An updated plotting_data tibble data frame
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' obj <- setPlottingData(obj, updated_plotting_data)
#' }

setGeneric("setPlottingData",
           function(MetAlyzer, plotting_data)
             standardGeneric("setPlottingData")
)
setMethod("setPlottingData",
          "MetAlyzer",
          function(MetAlyzer, plotting_data) {
            MetAlyzer@plotting_data <- plotting_data
            return(MetAlyzer)
          }
)
