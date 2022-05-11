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
             standardGeneric("metabolites")
)

#' @describeIn metabolites Get metabolites
setMethod("metabolites",
          "MetAlyzer",
          function(object) {
            object@metabolites[["filtered"]]
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
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' meta_data <- metaData(data)
#' head(meta_data)
#' }

setGeneric("metaData",
           function(object)
             standardGeneric("metaData")
)

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
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' raw_data <- rawData(data)
#' head(raw_data, c(5, 5))
#' }

setGeneric("rawData",
           function(object)
             standardGeneric("rawData")
)

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
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#'
#' quant_status <- quantStatus(data)
#' head(quant_status, c(5, 5))
#' }

setGeneric("quantStatus",
           function(object)
             standardGeneric("quantStatus")
)

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
          }
)
