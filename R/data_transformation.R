#' @title Transformation
#'
#' @description This function performs transformation of imputed concentration values
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
#' @param metalyzer_se a MetAlyzer object
#' @return An updated aggregated_data tibble data frame
#' @import dplyr
#' @importFrom rlang .data
#' 
#' @keywords internal
data_transformation <- function(metalyzer_se) {
  aggregated_data <- metalyzer_se@metadata$aggregated_data
  aggregated_data <- mutate(aggregated_data,
                            log2_Conc = transform(.data$imputed_Conc, base::log2),
                            .after = .data$imputed_Conc)
  metalyzer_se@metadata$aggregated_data <- aggregated_data
  return(metalyzer_se)
}
