#' Filter meta data
#'
#' This function updates the "Filter" column in meta_data to filter out samples.
#'
#' @param object MetAlyzer object
#' @param column A length-one character vector specifying which column of
#' meta_data to use for filtering
#' @param keep A vector specifying which samples to keep
#' @param remove A vector specifying which samples to remove
#'
#' @keywords internal

filter_meta_data <- function(object, column, keep, remove) {
  old_filter <- object@meta_data$Filter
  if (!is.null(keep)) {
    new_filter <- object@meta_data[,column] %in% keep
  } else if (!is.null(remove)) {
    new_filter <- ! object@meta_data[,column] %in% remove
  }
  updated_filter <- sapply(1:nrow(object@meta_data), function(i) {
    old_filter[i] & new_filter[i]
  })
  object@meta_data$Filter <- updated_filter
  return(object)
}
