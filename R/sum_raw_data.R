#' Summarize raw data
#'
#' This function prints quantiles and NAs of raw data.
#'
#' @param object MetAlyzer object
#'
#' @importFrom stats quantile
#'
#' @keywords internal

sum_raw_data <- function(object) {
  raw_data <- rawData(object)
  nas <- sum(is.na(raw_data))
  total <- nrow(raw_data) * ncol(raw_data)
  n_nas <- colSums(is.na(raw_data))
  na_metabolites <- colnames(raw_data)[n_nas > 0]

  cat("-------------------------------------\n")
  cat("Quantiles:\n")
  print(stats::quantile(raw_data, na.rm = TRUE))
  cat(paste0("\nNAs: ", nas, " (", round(nas/total*100, 2),"%)\n"))
  cat("-------------------------------------\n")

  return(na_metabolites)
}
