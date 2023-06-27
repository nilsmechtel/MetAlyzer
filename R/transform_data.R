#' Transformation
#'
#' This function performs transformation of imputed concentration values
#' (imputed_Conc).
#'
#' @param vec MetAlyzer metalyzer
#' @param func A function for transformation
#'
#' @keywords internal

transform <- function(vec, func) {
  vec[vec > 0 & !is.na(vec)] <- func(vec[vec > 0 & !is.na(vec)])
  return(vec)
}

#' @title Transform aggregated data
#'
#' @description This function performs a transformation of imputed concentration
#' values (imputed_Conc) with a given function. NA values are skipped. The
#' transformed values are added to aggregated_data in an extra column
#' transf_Conc.
#'
#' @param aggregated_data aggregated_data tibble data frame
#' @return An updated aggregated_data tibble data frame
#' @export
transform_data <- function(aggregated_data) {
  aggregated_data <- mutate(aggregated_data,
                            log2_Conc = transform(imputed_Conc, base::log2),
                            .after = imputed_Conc)
  return(aggregated_data)
}