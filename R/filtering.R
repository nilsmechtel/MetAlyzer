#' Filter metabolites
#'
#' This function filters out certain classes of the metabolites vector and raw_data and quant_status columns
#' @param object MetAlyzer object
#' @param class_name A class to be filtered out
#' 
#' @return 
#' @export
#' 
#' @examples

filter_classes <- function(object, class_name="Metabolism Indicators") {
  if (class_name %in% names(object@metabolites)) {
    not_indicators <- which(names(object@metabolites) != class_name)
    metabolites <- object@metabolites[not_indicators]
    object@metabolites <- metabolites
    print("-------------------------------------")
    print(paste(class_name, "were removed from metabolites."))
    if (nrow(object@raw_data) > 0) {
      object@raw_data <- object@raw_data[, metabolites]
      print(paste(class_name, "were removed from raw_data."))
    }
    if (nrow(object@quant_status) > 0) {
      object@quant_status <- object@quant_status[, metabolites]
      print(paste(class_name, "were removed from quant_status."))
    }
  } else {
    print("-------------------------------------")
    print(paste("No", class_name, "to filter! Returning original object."))
  }
  return(object)
}