#' Threshold CV
#'
#' This function assigns a CV value according to a vector of thresholds.
#' @param x A CV value
#' @param cv_threshs A numeric vector of thresholds
#'
#' @keywords internal
set_threshold <- function(x, cv_threshs) {
  levels <- c(paste0("max", cv_threshs*100),
              paste0("more", utils::tail(cv_threshs, 1)*100))
  cv_threshs <- sort(cv_threshs, decreasing = TRUE)
  if (is.na(x)) {
    v <- NA
  } else {
    for (t in cv_threshs) {
      if (x <= t) {
        v <- paste0("max", t*100)
      }
    }
    if (x > cv_threshs[1]) {
      v <- paste0("more", cv_threshs[1]*100)
    }
  }
  v <- factor(v, levels = levels)
  return(v)
}

#' @title Add mean, SD and CV
#'
#' @description This function calculates the mean, standard deviation (SD)
#' and the coefficient of variation (CV) for each group and adds them to
#' aggregated_data.
#'
#' @param aggregated_data aggregated_data tibble data frame
#' @param cv_thresholds A numeric vector of upper thresholds (CV <= t) between 0
#' and 1 for CV categorization.
#' @return An updated aggregated_data tibble data frame
#' @import dplyr
#' @export
calc_CV <- function(aggregated_data, cv_thresholds = c(0.1, 0.2, 0.3)) {
  grouping_vars <- groups(aggregated_data)
  cat(paste0("Calculate mean and coefficient of variation (groupwise: ",
             paste(grouping_vars, collapse = " * "), ")...  "))
  aggregated_data <- aggregated_data %>%
    dplyr::mutate(Mean = mean(Concentration, na.rm = TRUE),
           SD = stats::sd(Concentration, na.rm = TRUE),
           CV = SD / Mean,
           CV_thresh = sapply(CV, function(x) {
             set_threshold(x, cv_threshs = cv_thresholds)
           }),
           .after = Concentration)
  cat("finished!\n")
  return(aggregated_data)
}
