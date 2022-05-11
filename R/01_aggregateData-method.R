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
#' @return An updated MetAlyzer object
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
           function(object, ...,
                    id_col=NULL,
                    cv_threshs=c(0.1, 0.2, 0.3),
                    valid_vec=c("Valid", "LOQ"),
                    ko_vec=c(""),
                    valid_thresh=0.5,
                    impute_to=0.2,
                    impute_NA=FALSE)
             standardGeneric("aggregateData")
)

#' @describeIn aggregateData Aggregate data
setMethod("aggregateData",
          "MetAlyzer",
          function(object, ..., id_col, cv_threshs, valid_vec, valid_thresh,
                   impute_to, impute_NA) {
            aggregate_data(object, ...,
                           id_col=deparse(substitute(id_col)),
                           cv_threshs=cv_threshs,
                           valid_vec=valid_vec,
                           ko_vec=ko_vec,
                           valid_thresh=valid_thresh,
                           impute_to=impute_to,
                           impute_NA=impute_NA)
          }
)
