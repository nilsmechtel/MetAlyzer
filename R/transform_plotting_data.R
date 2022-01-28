#' Transform plotting data
#'
#' This method performs a transformation of imputed concentration values
#' (imp_Conc) with a given function. NA values are skipped. The transformed
#' values are added to plotting_data in an extra column transf_Conc.
#'
#' @param object MetAlyzer object
#' @param func A function to transform concentration values with
#'
#' @importFrom rlang .data
#'
#' @keywords internal

transform_plotting_data <- function(object, func) {
  object@plotting_data <- mutate(object@plotting_data,
                                 transf_Conc = transform(.data$imp_Conc, func),
                                 .after = .data$imp_Conc)
  return(object)
}


#' Transformation
#'
#' This function performs transformation of imputed concentration values
#' (imp_Conc).
#'
#' @param vec MetAlyzer object
#' @param func A function for transformation
#'
#' @keywords internal

transform <- function(vec, func) {
  vec[vec > 0 & !is.na(vec)] <- func(vec[vec > 0 & !is.na(vec)])
  return(vec)
}
