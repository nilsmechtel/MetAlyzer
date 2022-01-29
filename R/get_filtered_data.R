#' Get filtered data
#'
#' This function returns the filtered meta_data, raw_data or quant_status
#' data frame.
#'
#' @param object MetAlyzer object
#' @param slot A character value specifying which data frame to slice
#'
#' @import dplyr
#'
#' @keywords internal

get_filtered_data <- function(object, slot) {
  if (slot == "meta") {
    if (nrow(object@meta_data) > 0) {
      sliced_df <- object@meta_data %>%
        filter(Filter) %>%
        select(-Filter)
    } else {
      sliced_df <- object@meta_data
    }
  } else if (slot == "data") {
    if (nrow(object@raw_data) > 0) {
      sliced_df <- object@raw_data[object@meta_data$Filter,
                                   object@metabolites]
    } else {
      sliced_df <- object@raw_data
    }
  } else if (slot == "quant") {
    if (nrow(object@quant_status) > 0) {
      sliced_df <- object@quant_status[object@meta_data$Filter,
                                       object@metabolites]
    } else {
      sliced_df <- object@quant_status
    }
  }
  return(sliced_df)
}
