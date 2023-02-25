# MetAlyzer ---------------------------------------------------------------

#' A S4 class to read and analyze 'MetIDQ' output
#'
#' @slot file_path A length-one character vector giving the file path
#' @slot sheet A length-one numeric vector giving the sheet index
#' @slot metabolites A character vector with all 630 measured metabolites and
#' optional 234 additional metabolism indicators
#' @slot meta_data A data frame containing any meta data;
#' dimension: # samples x # meta variables
#' @slot raw_data A data frame containing all raw measurements;
#' dimension: # samples x # metabolites
#' @slot quant_status A data frame containing the quantification status of all
#' measurements; dimension: # samples x # metabolites
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzer_dataset(file_path = fpath)
#' }

MetAlyzer <- setClass(
  "MetAlyzer",
  slots = list(
    file_path = "character",
    sheet = "numeric",
    metabolites = "list",
    meta_data = "data.frame",
    raw_data = "data.frame",
    quant_status = "data.frame"
  )
)


# access data -------------------------------------------------------------

#' Get metabolites
#'
#' This method returns the filtered metabolites vector.
#'
#' @param object MetAlyzer object
#'
#' @return The metabolites vector
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' metabolites(data)
#' }

setGeneric("metabolites",
           function(object)
             standardGeneric("metabolites"))

#' @describeIn metabolites Get metabolites
setMethod("metabolites",
          "MetAlyzer",
          function(object) {
            object@metabolites[["filtered"]]
          })


#' Get meta data
#'
#' This method returns the meta_data data frame with filtered samples (rows) and
#' all columns.
#'
#' @param object MetAlyzer object
#'
#' @return The meta_data data frame
#'
#' @importFrom dplyr filter select
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' meta_data <- metaData(data)
#' head(meta_data)
#' }

setGeneric("metaData",
           function(object)
             standardGeneric("metaData"))

#' @describeIn metaData Get meta data
setMethod("metaData",
          "MetAlyzer",
          function(object) {
            if (nrow(object@meta_data) > 0) {
              meta_data <- object@meta_data %>%
                filter(Filter) %>%
                select(-Filter)
            } else {
              meta_data <- object@meta_data
            }
            return(meta_data)
          })


#' Get raw data
#'
#' This method returns the raw_data data frame with filtered samples (rows) and
#' metabolites (columns).
#'
#' @param object MetAlyzer object
#'
#' @return The raw_data data frame
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' raw_data <- rawData(data)
#' head(raw_data, c(5, 5))
#' }

setGeneric("rawData",
           function(object)
             standardGeneric("rawData"))

#' @describeIn rawData Get raw data
setMethod("rawData",
          "MetAlyzer",
          function(object) {
            if (nrow(object@raw_data) > 0) {
              raw_data <- object@raw_data[object@meta_data$Filter,
                                          metabolites(object)]
            } else {
              raw_data <- object@raw_data
            }
            return(raw_data)
          })


#' Get quantification status
#'
#' This method returns the quant_status data frame with filtered samples (rows)
#' and metabolites (columns).
#'
#' @param object MetAlyzer object
#'
#' @return The quant_status data frame
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' quant_status <- quantStatus(data)
#' head(quant_status, c(5, 5))
#' }

setGeneric("quantStatus",
           function(object)
             standardGeneric("quantStatus"))

#' @describeIn quantStatus Get quantification status
setMethod("quantStatus",
          "MetAlyzer",
          function(object) {
            if (nrow(object@quant_status) > 0) {
              quant_status <- object@quant_status[object@meta_data$Filter,
                                                  metabolites(object)]
            } else {
              quant_status <- object@quant_status
            }
            return(quant_status)
          })


# display data ------------------------------------------------------------

#' Show a 'MetAlyzer' object
#'
#' This method shows a summary of 'MetAlyzer' slot values.
#'
#' @param object MetAlyzer object
#'
#' @return A summary of the MetAlyzer object
#'
#' @importFrom methods show
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' show(data)
#' # or
#' data
#' }

setMethod("show",
          "MetAlyzer",
          function(object) {
            show_obj(object)
          })


#' Summarize raw data
#'
#' This method lists the number of each quantification status and its
#' percentage.
#'
#' @param object MetAlyzer object
#'
#' @return A summary of the quantification status
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' summarizeRawData(data)
#' }

setGeneric("summarizeRawData",
           function(object)
             standardGeneric("summarizeRawData"))

#' @describeIn summarizeRawData Summarize raw data
setMethod("summarizeRawData",
          "MetAlyzer",
          function(object) {
            sum_raw_data(object)
          })

#' Summarize quantification status
#'
#' This method lists the number of each quantification status and its
#' percentage.
#'
#' @param object MetAlyzer object
#'
#' @return A summary of the quantification status
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' summarizeQuantData(data)
#' }

setGeneric("summarizeQuantData",
           function(object)
             standardGeneric("summarizeQuantData"))

#' @describeIn summarizeQuantData Summarize quantification status
setMethod("summarizeQuantData",
          "MetAlyzer",
          function(object) {
            sum_quant_data(object)
          })


# manage meta data --------------------------------------------------------

#' Filter meta data
#'
#' This function updates the "Filter" column in meta_data to filter out samples.
#'
#' If both "keep" and "remove" arguments are used "keep" overwrites the
#' "remove" argument.
#'
#' @param object MetAlyzer object
#' @param column A column of meta_data for filtering
#' @param keep A vector defining which entries to keep from meta_data
#' @param remove A vector defining which entries to remove meta_data
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' data <- filterMetaData(data, column = Group, keep = 1:6)
#' # or
#' data <- filterMetaData(data, column = Group, remove = 7)
#' }

setGeneric("filterMetaData",
           function(object,
                    column,
                    keep = NULL,
                    remove = NULL)
             standardGeneric("filterMetaData"))

#' @describeIn filterMetaData Filter meta data
setMethod("filterMetaData",
          "MetAlyzer",
          function(object, column, keep, remove) {
            filter_meta_data(object, deparse(substitute(column)), keep, remove)
          })

#' Reset meta data
#'
#' This method resets the filter of meta_data.
#'
#' @param object MetAlyzer object
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' data <- resetMetaData(data)
#' }

setGeneric("resetMetaData",
           function(object)
             standardGeneric("resetMetaData"))

#' @describeIn resetMetaData Reset meta data
setMethod("resetMetaData",
          "MetAlyzer",
          function(object) {
            object@meta_data$Filter <- TRUE
            return(object)
          })

#' Update meta data
#'
#' This method adds another column to filtered meta_data.
#'
#' @param object MetAlyzer object
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
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' data <- updateMetaData(data, name = Date, new_colum = Sys.Date())
#' data <- updateMetaData(data, name = Analyzed, new_colum = TRUE)
#' }

setGeneric("updateMetaData",
           function(object, name, new_colum)
             standardGeneric("updateMetaData"))

#' @describeIn updateMetaData Update meta data
setMethod("updateMetaData",
          "MetAlyzer",
          function(object, name, new_colum) {
            update_meta_data(object, deparse(substitute(name)), new_colum)
          })

#' Rename meta data
#'
#' This method renames a column of meta_data using rename 'dplyr'.
#'
#' @param object MetAlyzer object
#' @param ... Use new_name = old_name to rename selected variables
#'
#' @return An updated MetAlyzer object
#'
#' @importFrom dplyr rename
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' data <- renameMetaData(data, Method = Group)
#' }

setGeneric("renameMetaData",
           function(object, ...)
             standardGeneric("renameMetaData"))

#' @describeIn renameMetaData Rename meta data
setMethod("renameMetaData",
          "MetAlyzer",
          function(object, ...) {
            object@meta_data <- rename(object@meta_data, ...)
            return(object)
          })


# manage metabolites ------------------------------------------------------

#' Filter metabolites
#'
#' This method filters out certain classes or metabolites of the metabolites
#' vector.
#' Note: "metabo_vec" overwrites the "class_name" argument!
#'
#' @param object MetAlyzer object
#' @param class_name A character value defining the class to be removed
#' @param metabo_vec A character vector defining metabolites to be removed
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' data <- filterMetabolites(data, class_name = "Metabolism Indicators")
#' # or
#' data <- filterMetabolites(data, metabo_vec = c("C0", "C2", "C3"))
#' }

setGeneric("filterMetabolites",
           function(object,
                    class_name = "Metabolism Indicators",
                    metabo_vec = NULL)
             standardGeneric("filterMetabolites"))

#' @describeIn filterMetabolites Filter metabolites
setMethod("filterMetabolites",
          "MetAlyzer",
          function(object, class_name, metabo_vec) {
            filter_metabolites(object, class_name, metabo_vec)
          })


#' Reset metabolites
#'
#' This method resets the filtering of metabolites.
#'
#' @param object MetAlyzer object
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' data <- resetMetabolites(data)
#' }

setGeneric("resetMetabolites",
           function(object)
             standardGeneric("resetMetabolites"))

#' @describeIn resetMetabolites Reset metabolites
setMethod("resetMetabolites",
          "MetAlyzer",
          function(object) {
            reset_metabolites(object)
          })


#' Drop invalid metabolites
#'
#' This method extracts metabolites with a given percentage of valid replicates
#' from aggregated_data and filters metabolites.
#'
#' @param object MetAlyzer object
#' @param aggregated_data aggregated_data tibble data frame
#' @param filter_col A boolean column from aggregated_data defines the filter
#' @param replicate_t A numeric lower threshold between 0 and 1 (t < x) to
#' determine valid metabolites with a given percentage of valid replicates.
#'
#' @return An updated MetAlyzer object
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' # To see an example, please check out the vignette.


setGeneric("dropInvalidMetabolites",
           function(object,
                    aggregated_data,
                    filter_col = .data$Valid_Replicates,
                    replicate_t = 0)
             standardGeneric("dropInvalidMetabolites"))

#' @describeIn dropInvalidMetabolites Drop invalid metabolites
setMethod("dropInvalidMetabolites",
          "MetAlyzer",
          function(object,
                   aggregated_data,
                   filter_col,
                   replicate_t) {
            invalid_metabolites <- get_invalid_metabolites(aggregated_data,
                                                           deparse(substitute(filter_col)),
                                                           replicate_t)
            object <-
              filter_metabolites(object, metabo_vec = invalid_metabolites)
            return(object)
          })


# export data -------------------------------------------------------------

#' Export filtered raw data as csv
#'
#' This method exports the filtered raw data in the CSV format.
#'
#' @param object MetAlyzer object
#' @param sample_id A column from meta_data for sample identification
#' @param file_path file path
#'
#' @export
#'
#' @examples
#' \dontrun{
#' print(1)
#' }

setGeneric("exportRawData",
           function(object, sample_id, file_path)
             standardGeneric("exportRawData"))

#' @describeIn exportRawData Export filtered raw data as csv
setMethod("exportRawData",
          "MetAlyzer",
          function(object, sample_id, file_path) {
            export_raw_data(object, deparse(substitute(sample_id)), file_path)
          })


# aggregate data ----------------------------------------------------------

#' Aggregate data
#'
#' This method reshapes raw_data, quant_status and meta_data and combines them
#' together with basic statistics in a tibble data frame. "aggregated_data" is
#' grouped by metabolites as well as the selection of additional variables.
#' Statistics are then calculated for each group. The format is suitable for
#' plotting with ggplot2.
#'
#' @param object MetAlyzer object
#' @param ... A selection of columns from meta_data to add to aggregated data
#' frame that will be used as extra grouping variables.
#' @param id_col Sample ID column from meta_data that is added to aggregated data
#' frame but is not used as grouping variable.
#' @param cv_threshs A numeric vector of upper thresholds (CV <= t) between 0
#' and 1 for CV categorization.
#' @param valid_vec A character vector containing each quantification status that
#' is considered to be a valid measurement.
#' @param ko_vec A character vector containing each quantification status that
#' invalidates a set of repetitions once it is encountered.
#' @param valid_thresh A numeric lower threshold between 0 and 1 (t < x) to
#' determine valid replicates based on their consideration of a valid measurement.
#' @param impute_to A numeric value below 1
#' @param impute_NA A logical value whether to impute NA values
#'
#' @return A tibble data frame
#'
#' @importFrom rlang .data
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzer_dataset(file_path = fpath)
#'
#' data <- aggregateData(data, Tissue, Group,
#'                       ungrouped = NULL,
#'                       ts = c(0.1, 0.2, 0.3),
#'                       valid_vec = c("Valid", "LOQ"), t = 0.5)
#' }

setGeneric("aggregateData",
           function(object,
                    ...,
                    id_col = NULL,
                    cv_threshs = c(0.1, 0.2, 0.3),
                    valid_vec = c("Valid", "LOQ"),
                    ko_vec = c(""),
                    valid_thresh = 0.5,
                    impute_to = 0.2,
                    impute_NA = FALSE)
             standardGeneric("aggregateData"))

#' @describeIn aggregateData Aggregate data
setMethod("aggregateData",
          "MetAlyzer",
          function(object,
                   ...,
                   id_col,
                   cv_threshs,
                   valid_vec,
                   ko_vec,
                   valid_thresh,
                   impute_to,
                   impute_NA) {
            aggregate_data(
              object,
              ...,
              id_col = deparse(substitute(id_col)),
              cv_threshs = cv_threshs,
              valid_vec = valid_vec,
              ko_vec = ko_vec,
              valid_thresh = valid_thresh,
              impute_to = impute_to,
              impute_NA = impute_NA
            )
          })
