#' Impute plotting data
#'
#' This function imputes zero concentration values (Concentration) with the
#' minimal positive value multiplied by i. If all values are zero or NA, they
#' are set to NA. The imputed values are added to plotting_data in an extra
#' column imp_Conc.
#'
#' @param object MetAlyzer object
#' @param ... Variables to group by
#' @param i A numeric value below 1
#' @param imputeNA Logical value whether to impute NA values
#'
#' @importFrom rlang .data
#'
#' @keywords internal

impute_plotting_data <- function(object, i, imputeNA, ...) {
  plotting_data <- object@plotting_data
  grouping_vars <- group_vars(plotting_data)
  plotting_data <- plotting_data %>%
    group_by(...) %>%
    mutate(imp_Conc = zero_imputation(.data$Concentration, i, imputeNA),
           .after = .data$Concentration) %>%
    group_by_at(grouping_vars)
  object@plotting_data <- plotting_data
  return(object)
}


#' Zero imputation
#'
#' This function performs zero imputation with the minimal positive value times
#' i.
#'
#' @param vec A numeric vector containing the concentration values
#' @param i A numeric value below 1)
#' @param imputeNA Logical value whether to impute NA values
#'
#' @keywords internal

zero_imputation <- function(vec, i, imputeNA) {
  non_zero <- vec[vec > 0 & !is.na(vec)]
  imp_v <- ifelse(length(non_zero) > 0, min(non_zero) * i, NA)
  vec[vec == 0] <- imp_v
  if (imputeNA) {
    vec[is.na(vec)] <- imp_v
  }
  return(vec)
}
