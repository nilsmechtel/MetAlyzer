#' Filter metabolites
#'
#' This function filters out certain classes or metabolites of the metabolites
#' vector
#' @param object MetAlyzer object
#' @param class_name A class to be filtered out
#' @param metabo_vec A character vector defining metabolites to be removed
#'
#' @keywords internal

filter_metabolites <- function(object, class_name, metabo_vec) {
  metabolites <- object@metabolites
  orig_length <- length(metabolites)
  if (!is.null(metabo_vec)) { # if metabo_vec is given
    if (any(metabo_vec %in% metabolites)) {
      metabolites <- metabolites[which(!metabolites %in% metabo_vec)]
    }
  } else if (class_name %in% names(metabolites)) {
    metabolites <- metabolites[which(names(metabolites) != class_name)]
  }
  diff <- orig_length - length(metabolites)
  if (diff == 1) {
    cat("1 metabolite was filtered!\n")
  } else if (diff > 1) {
    cat(paste(diff, "metabolites were filtered!\n"))
  } else {
    cat("No metabolites were filtered!\n")
  }
  object@metabolites <- metabolites
  return(object)
}


#' Filter meta data
#'
#' This function updates the "Filter" column in meta_data to filter out samples
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

#' Update meta data
#'
#' This function adds another column to unfiltered meta_data
#' @param MetAlyzer MetAlyzer object
#' @param name The new column name
#' @param new_colum A vector for the new column (length has to be same as the
#' number of filtered samples)
#'
#' @keywords internal

update_meta_data <- function(MetAlyzer, name, new_colum) {
  if (class(new_colum) == "factor") {
    levels <- levels(new_colum)
  } else {
    levels <- unique(new_colum)
  }
  meta_data <- MetAlyzer@meta_data
  chr_name <- deparse(substitute(name))
  meta_data[,chr_name] <- factor(NA, levels = levels)
  meta_data[,chr_name][meta_data$Filter == TRUE] <- new_colum
  MetAlyzer@meta_data <- meta_data
  return(MetAlyzer)
}


#' Get filtered data
#'
#' This function filters meta_data, raw_data or quant_status
#' @param object MetAlyzer object
#' @param slot A length-one character vector specifying which data frame to slice
#' @param verbose If TRUE prints which filtered data frame is returned
#'
#' @import dplyr
#'
#' @keywords internal

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
