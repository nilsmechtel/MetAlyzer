#' Filter metabolites
#'
#' This function filters out certain classes or metabolites of the metabolites
#' vector.
#'
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
