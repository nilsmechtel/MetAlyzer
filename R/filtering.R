#' Filter metabolites
#'
#' This function filters out certain classes of the metabolites vector
#' @param object MetAlyzer object
#' @param class_name A class to be filtered out
#'
#' @return
#'
#' @export
#'
#' @examples

filter_metabolites <- function(object, class_name="Metabolism Indicators") {
  cat("-------------------------------------\n")
  if (class_name %in% names(object@metabolites)) {
    new_metabolites <- which(names(object@metabolites) != class_name)
    object@metabolites <- object@metabolites[new_metabolites]
    if (! class_name %in% names(object@metabolites)) {
      cat(paste(class_name, "were removed from metabolites.\n"))
    }
  } else {
    cat(paste("No", class_name, "to filter!\n"))
  }
  return(object)
}


#' Filter meta data
#'
#' This function updates the "Filter" column in meta_data to filter out samples
#' @param object MetAlyzer object
#' @param column A length-one character vector specifying which column of
#' meta_data to use for filtering
#' @param keep A character vector specifying which samples to keep
#' @param remove A character vector specifying which samples to remove
#'
#' @return
#'
#' @export
#'
#' @examples

filter_meta_data <- function(object, column, keep=NULL, remove=NULL) {
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


#' Get filtered data
#'
#' This function filters meta_data, raw_data or quant_status
#' @param object MetAlyzer object
#' @param slot A length-one character vector specifying which data frame to slice
#'
#' @return
#'
#' @import dplyr
#' @export
#'
#' @examples

get_filtered_data <- function(object, slot, verbose=TRUE) {
  if (slot == "meta") {
    if (nrow(object@meta_data) > 0) {
      sliced_df <- object@meta_data %>%
        filter(Filter) %>%
        select(-Filter)
      if (verbose) {
        cat("-------------------------------------\n")
        cat("Returning filtered meta data\n")
      }
    } else {
      sliced_df <- object@meta_data
    }
  } else if (slot == "data") {
    if (nrow(object@raw_data) > 0) {
      sliced_df <- object@raw_data[object@meta_data$Filter, object@metabolites]
      if (verbose) {
        cat("-------------------------------------\n")
        cat("Returning filtered raw data\n")
      }
    } else {
      sliced_df <- object@raw_data
    }
  } else if (slot == "quant") {
    if (nrow(object@quant_status) > 0) {
      sliced_df <- object@quant_status[object@meta_data$Filter, object@metabolites]
      if (verbose) {
        cat("-------------------------------------\n")
        cat("Returning filtered quantification status\n")
      }
    } else {
      sliced_df <- object@quant_status
    }
  }
  return(sliced_df)
}
