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
             standardGeneric("exportRawData")
)

#' @describeIn exportRawData Export filtered raw data as csv
setMethod("exportRawData",
          "MetAlyzer",
          function(object, sample_id, file_path) {
            export_raw_data(object, deparse(substitute(sample_id)), file_path)
          }
)
