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
                    class_name="Metabolism Indicators",
                    metabo_vec=NULL)
             standardGeneric("filterMetabolites")
)

#' @describeIn filterMetabolites Filter metabolites
setMethod("filterMetabolites",
          "MetAlyzer",
          function(object, class_name, metabo_vec) {
            filter_metabolites(object, class_name, metabo_vec)
          }
)


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
             standardGeneric("resetMetabolites")
)

#' @describeIn resetMetabolites Reset metabolites
setMethod("resetMetabolites",
          "MetAlyzer",
          function(object) {
            reset_metabolites(object)
          }
)


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
                    filter_col=.data$Valid_Replicates,
                    replicate_t=0)
             standardGeneric("dropInvalidMetabolites")
)

#' @describeIn dropInvalidMetabolites Drop invalid metabolites
setMethod("dropInvalidMetabolites",
          "MetAlyzer",
          function(object, aggregated_data, filter_col, replicate_t) {
            invalid_metabolites <- get_invalid_metabolites(aggregated_data,
                                                           deparse(substitute(filter_col)),
                                                           replicate_t)
            object <- filter_metabolites(object, metabo_vec=invalid_metabolites)
            return(object)
          }
)
