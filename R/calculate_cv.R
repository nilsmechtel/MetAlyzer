#' @title Add mean, SD and CV
#'
#' @description This function calculates the mean, standard deviation (SD)
#' and the coefficient of variation (CV) for each group and adds them to
#' aggregated_data.
#'
#' @param metalyzer_se A Metalyzer object
#' @param groups A vector of column names of aggregated_data to calculate mean,
#' SD and CV group wise. If the column does not exists in aggregated_data it is
#' automatically added from meta data. The default value is set to NULL, which
#' uses the existing grouping of aggregated_data. 
#' @param cv_thresholds A numeric vector of upper thresholds (CV <= t) between 0
#' and 1 for CV categorization.
#' @param na.rm a logical evaluating to TRUE or FALSE indicating whether NA
#' values should be stripped before the computation proceeds.
#' @return An updated aggregated_data tibble data frame
#' 
#' @import dplyr
#' @importFrom rlang .data
#' @export
#' 
#' @examples
#' metalyzer_se <- MetAlyzer_dataset(file_path = example_extraction_data())
#' metalyzer_se <- renameMetaData(
#'   metalyzer_se,
#'   Extraction_Method = "Sample Description"
#' )
#' metalyzer_se <- filterMetaData(
#'   metalyzer_se,
#'   Tissue == "Drosophila"
#' )
#' metalyzer_se <- calculate_cv(
#'   metalyzer_se,
#'   groups = c("Tissue", "Extraction_Method", "Metabolite"),
#'   cv_thresholds = c(0.1, 0.2, 0.3),
#'   na.rm = TRUE
#' )
calculate_cv <- function(metalyzer_se,
                         groups = NULL,
                         cv_thresholds = c(0.1, 0.2, 0.3),
                         na.rm = TRUE) {
  if (is.null(groups)) {
    groups <- as.character(groups(metalyzer_se@metadata$aggregated_data))
  } else {
    for (group in groups) {
      if (!group %in% colnames(metalyzer_se@metadata$aggregated_data)) {
        metalyzer_se <- expand_aggregated_data(metalyzer_se, meta_data_column = group)
      }
    }
  }
  aggregated_data <- metalyzer_se@metadata$aggregated_data %>%
    group_by_at(unique(c(groups, "Metabolite")))
  
  cat(paste0("Info: Calculating mean and coefficient of variation (groupwise: ",
             paste(groups(aggregated_data), collapse = " * "), ")...  "))
  aggregated_data <- aggregated_data %>%
    dplyr::mutate(Mean = mean(.data$Concentration, na.rm = na.rm),
           SD = stats::sd(.data$Concentration, na.rm = na.rm),
           CV = .data$SD / .data$Mean,
           CV_thresh = sapply(.data$CV, function(x) {
             set_threshold(x, cv_threshs = cv_thresholds)
           }),
           .after = .data$Concentration)
  aggregated_data$CV[is.na(aggregated_data$CV)] <- NA
  cat("finished!\n")
  metalyzer_se@metadata$aggregated_data <- aggregated_data
  return(metalyzer_se)
}


#' @title Threshold CV
#'
#' @description This function assigns a CV value according to a vector of thresholds.
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
