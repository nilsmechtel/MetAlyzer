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
          }
)


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
             standardGeneric("summarizeRawData")
)

#' @describeIn summarizeRawData Summarize raw data
setMethod("summarizeRawData",
          "MetAlyzer",
          function(object) {
            sum_raw_data(object)
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
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' summarizeQuantData(data)
#' }

setGeneric("summarizeQuantData",
           function(object)
             standardGeneric("summarizeQuantData")
)

#' @describeIn summarizeQuantData Summarize quantification status
setMethod("summarizeQuantData",
          "MetAlyzer",
          function(object) {
            sum_quant_data(object)
          }
)
