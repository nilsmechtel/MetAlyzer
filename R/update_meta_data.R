#' Update meta data
#'
#' This function adds another column to filtered meta_data.
#'
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
  meta_data[,name] <- factor(NA, levels = levels)
  meta_data[,name][meta_data$Filter == TRUE] <- new_colum
  MetAlyzer@meta_data <- meta_data
  return(MetAlyzer)
}
