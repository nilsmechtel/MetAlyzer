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
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' show(obj)
#' # or
#' obj


setMethod("show",
          "MetAlyzer",
          function(object) {
            show_obj(object)
          }
)


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
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' summariseQuantData(obj)

setGeneric("summariseQuantData",
           function(object)
             standardGeneric("summariseQuantData")
)

#' @describeIn summariseQuantData Summarize quantification status
setMethod("summariseQuantData",
          "MetAlyzer",
          function(object) {
            sum_quant_data(object)
          }
)


#' Get meta data
#'
#' This method returns the meta_data data frame with filtered samples (rows) and
#' all columns.
#'
#' @param object MetAlyzer object
#'
#' @return The meta_data data frame
#'
#' @export
#'
#' @examples
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' meta_data <- metaData(obj)
#' head(meta_data)

setGeneric("metaData",
           function(object)
             standardGeneric("metaData")
)

#' @describeIn metaData Get meta data
setMethod("metaData",
          "MetAlyzer",
          function(object) {
            get_filtered_data(object, "meta")
          }
)


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
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' metabolites(obj)

setGeneric("metabolites",
           function(object)
             standardGeneric("metabolites")
)

#' @describeIn metabolites Get metabolites
setMethod("metabolites",
          "MetAlyzer",
          function(object) {
            object@metabolites
          }
)


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
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' raw_data <- rawData(obj)
#' head(raw_data, c(5, 5))

setGeneric("rawData",
           function(object)
             standardGeneric("rawData")
)

#' @describeIn rawData Get raw data
setMethod("rawData",
          "MetAlyzer",
          function(object) {
            get_filtered_data(object, "data")
          }
)


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
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' quant_status <- quantStatus(obj)
#' head(quant_status, c(5, 5))

setGeneric("quantStatus",
           function(object)
             standardGeneric("quantStatus")
)

#' @describeIn quantStatus Get quantification status
setMethod("quantStatus",
          "MetAlyzer",
          function(object) {
            get_filtered_data(object, "quant")
          }
)


#' Get plotting data
#'
#' This method returns the plotting_data tibble data frame.
#'
#' @param object MetAlyzer object
#'
#' @return The plotting_data data frame
#'
#' @export
#'
#' @examples
#' \donttest{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#' obj <- createPlottingData(obj, Tissue, Group)
#'
#' plottingData(obj)
#' }

setGeneric("plottingData",
           function(object)
             standardGeneric("plottingData")
)

#' @describeIn plottingData Get plotting data
setMethod("plottingData",
          "MetAlyzer",
          function(object) {
            object@plotting_data
          }
)
