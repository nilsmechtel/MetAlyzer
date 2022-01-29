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
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' obj <- filterMetaData(obj, column = Group, keep = 1:6)
#' # or
#' obj <- filterMetaData(obj, column = Group, remove = 7)

setGeneric("filterMetaData",
           function(object, column,
                    keep=NULL,
                    remove=NULL)
             standardGeneric("filterMetaData")
)

#' @describeIn filterMetaData Filter meta data
setMethod("filterMetaData",
          "MetAlyzer",
          function(object, column, keep, remove) {
            filter_meta_data(object, deparse(substitute(column)), keep, remove)
          }
)

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
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' obj <- resetMetaData(obj)

setGeneric("resetMetaData",
           function(object)
             standardGeneric("resetMetaData")
)

#' @describeIn resetMetaData Reset meta data
setMethod("resetMetaData",
          "MetAlyzer",
          function(object) {
            object@meta_data$Filter <- TRUE
            return(object)
          }
)

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
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' obj <- updateMetaData(obj, name = Date, new_colum = Sys.Date())
#' obj <- updateMetaData(obj, name = Analyzed, new_colum = TRUE)

setGeneric("updateMetaData",
           function(object, name, new_colum)
             standardGeneric("updateMetaData")
)

#' @describeIn updateMetaData Update meta data
setMethod("updateMetaData",
          "MetAlyzer",
          function(object, name, new_colum) {
            update_meta_data(object, deparse(substitute(name)), new_colum)
          }
)

#' Rename meta data
#'
#' This method renames a column of meta_data using rename 'dplyr'.
#'
#' @param object MetAlyzer object
#' @param ... Use new_name = old_name to rename selected variables
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' obj <- renameMetaData(obj, Method = Group)

setGeneric("renameMetaData",
           function(object, ...)
             standardGeneric("renameMetaData")
)

#' @describeIn renameMetaData Rename meta data
setMethod("renameMetaData",
          "MetAlyzer",
          function(object, ...) {
            object@meta_data <- rename(object@meta_data, ...)
            return(object)
          }
)
