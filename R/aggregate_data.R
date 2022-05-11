#' Aggregate_data
#'
#' This function reshapes raw_data, quant_status and meta_data and combines them
#' together with basic statistics in a tibble data frame. "aggregated_data" is
#' grouped by metabolites as well as the selection of additional variables.
#' Statistics are then calculated for each group. The format is suitable for
#' plotting with ggplot2.
#'
#' @param object MetAlyzer object
#' @param ... A selection of columns from meta_data to add to aggregated data
#' frame that will be used as extra grouping variables.
#' @param id_col Sample ID column from meta_data that is added to aggregated data
#' frame but is not used as grouping variable.
#' @param cv_threshs A numeric vector of upper thresholds (CV <= t) between 0
#' and 1 for CV categorization.
#' @param valid_vec A character vector containing each quantification status that
#' is considered to be a valid measurement.
#' @param ko_vec A character vector containing each quantification status that
#' invalidates a set of repetitions once it is encountered.
#' @param valid_thresh A numeric lower threshold between 0 and 1 (t < x) to
#' determine valid replicates based on their consideration of a valid measurement.
#' @param impute_to A numeric value below 1
#' @param impute_NA A logical value whether to impute NA values
#'
#' @keywords internal

aggregate_data <- function(object, ..., id_col, cv_threshs, valid_vec, ko_vec,
                           valid_thresh, impute_to, impute_NA) {
  cat("Reshape and merge data...  ")
  aggregated_data <- reshape_data(object, ..., id_col = id_col)
  cat("finished!\n")
  grouping_vars <- groups(aggregated_data)

  cat(paste0("Calculate mean and coefficient of variation (groupwise: ",
             paste(grouping_vars, collapse = " * "), ")...  "))
  aggregated_data <- calc_CV(aggregated_data, cv_threshs = cv_threshs)
  cat("finished!\n")
  cat(paste0("Calculate valid replicates (groupwise: ",
             paste(grouping_vars, collapse = " * "), ")...  "))
  aggregated_data <- valid_measurement(aggregated_data, valid_vec, ko_vec,
                                       valid_thresh)
  cat("finished!\n")
  cat("Impute concentrations (groupwise: Metabolite) with",
      paste0(round(impute_to * 100), "%"), "of the minimal positive value...  ")

  aggregated_data <- impute_data(aggregated_data, impute_to, impute_NA)
  cat("finished!\n")
  cat("Log2 transform imputed concentrations...  ")
  aggregated_data <- transform_plotting_data(aggregated_data)
  cat("finished!\n")
  return(aggregated_data)
}


#' Reshape data
#'
#' This function reshapes raw_data, quant_status and meta_data and combines them
#' in a tibble data frame for plotting with 'ggplot2'.
#' @param object MetAlyzer object
#' @param ... A selection of columns from meta_data to add to aggregated data frame
#' @param id_col Sample ID column from meta_data that is added to aggregated data
#' frame but is not used as grouping variable.
#'
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom rlang .data
#'
#' @keywords internal

reshape_data <- function(object, ..., id_col) {
  meta_data <- metaData(object)
  metabolites <- object@metabolites[["filtered"]]
  raw_data <- rawData(object)
  quant_status <- quantStatus(object)
  if (nchar(id_col) == 0 | id_col == "NULL") {
    id_col <- NULL
  }
  if (nrow(meta_data) > 0) {
    extra_columns <- select(meta_data, c(id_col, ...))
    raw_data <- bind_cols(raw_data, extra_columns)
    gathered_data <- gather(raw_data, key = "Metabolite",
                            value = "Concentration", -colnames(extra_columns))
    gathered_status <- gather(quant_status, key = "Metabolite",
                              value = "Status")
    col_names <- colnames(select(gathered_data, -id_col))
    group_cols <- col_names[1:which(col_names == "Metabolite")]
    aggregated_data <- gathered_data %>%
      group_by_at(group_cols) %>%
      mutate(Class = sapply(.data$Metabolite, function(x) {
        names(metabolites[metabolites == x])
      }),
      .after = .data$Metabolite)
    aggregated_data$Metabolite <- factor(aggregated_data$Metabolite,
                                         levels = unique(metabolites))
    aggregated_data$Class <- factor(aggregated_data$Class,
                                    levels = unique(names(metabolites)))
    aggregated_data$Status <- factor(gathered_status$Status,
                                     levels = levels(quant_status[,1]))
    for (col_name in colnames(extra_columns)) {
      col <- extra_columns[, col_name]
      if (inherits(col, "factor")) {
        levels <- levels(col)
      } else if (!any(grepl("\\D", col))) {
        levels <- unique(as.character(sort(as.numeric(col))))
      } else {
        levels <- unique(col)
      }
      aggregated_data[, col_name] <- factor(unlist(aggregated_data[,col_name]),
                                            levels = levels)
    }
    aggregated_data <- arrange_at(aggregated_data, c(group_cols, id_col))

    return(aggregated_data)
  } else {
    return(data.frame())
  }
}


#' Threshold CV
#'
#' This function assigns a CV value according to a vector of thresholds.
#' @param x A CV value
#' @param cv_threshs A numeric vector of thresholds
#'
#' @importFrom utils tail
#'
#' @keywords internal

set_threshold <- function(x, cv_threshs) {
  levels <- c(paste0("max", cv_threshs*100),
              paste0("more", tail(cv_threshs, 1)*100))
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


#' Add CV
#'
#' This function adds the mean, standard deviation (SD) and the
#' coefficient of variation (CV) to aggregated_data
#' @param aggregated_data aggregated_data tibble data frame
#' @param cv_threshs A numeric vector of thresholds
#'
#' @import dplyr
#' @importFrom stats sd
#' @importFrom rlang .data
#'
#' @keywords internal

calc_CV <- function(aggregated_data, cv_threshs) {
  aggregated_data <- aggregated_data %>%
    mutate(Mean = mean(.data$Concentration, na.rm = TRUE),
           SD = sd(.data$Concentration, na.rm = TRUE),
           CV = .data$SD / .data$Mean,
           CV_thresh = sapply(.data$CV, function(x) {
             set_threshold(x, cv_threshs = cv_threshs)
             }),
           .after = .data$Concentration)
  return(aggregated_data)
}


#' Add valid filter
#'
#' This function adds a filter column based on the quantification status. The
#' filter is True if the percentage of measurements with a quantification status
#' part of valid_vec is greater than the threshold valid_thresh.
#' @param aggregated_data aggregated_data tibble data frame
#' @param valid_vec A character vector containing each quantification status that
#' is considered to be a valid measurement
#' @param ko_vec A character vector containing each quantification status that
#' invalidates a set of repetitions once it is encountered.
#' @param valid_thresh A numeric threshold
#'
#' @import dplyr
#' @importFrom rlang .data
#'
#' @keywords internal

valid_measurement <- function(aggregated_data, valid_vec, ko_vec, valid_thresh) {
  filter_status <- function(vec) {
    status <- sum(vec %in% valid_vec) / length(vec) > valid_thresh
    if (status) {
      if (length(ko_vec) == 0) {
        status <- !any(vec %in% ko_vec)
      }
    }
    return(status)
  }
  aggregated_data <- mutate(aggregated_data,
                            Valid_Replicates = filter_status(.data$Status),
                            .after = .data$Status)
  return(aggregated_data)
}


#' Zero imputation
#'
#' This function performs zero imputation with the minimal positive value times
#' impute_to.
#'
#' @param vec A numeric vector containing the concentration values
#' @param impute_to A numeric value below 1
#' @param impute_NA Logical value whether to impute NA values
#'
#' @keywords internal

zero_imputation <- function(vec, impute_to, impute_NA) {
  non_zero <- vec[vec > 0 & !is.na(vec)]
  imp_v <- ifelse(length(non_zero) > 0, min(non_zero) * impute_to, NA)
  vec[vec == 0] <- imp_v
  if (impute_NA) {
    vec[is.na(vec)] <- imp_v
  }
  return(vec)
}


#' Impute aggregated data
#'
#' This function imputes zero concentration values (Concentration) with the
#' minimal positive value multiplied by impute_to. If all values are zero or NA,
#' they are set to NA. The imputed values are added to plotting_data in an extra
#' column imputed_Conc.
#'
#' @param aggregated_data aggregated_data tibble data frame
#' @param impute_to A numeric value below 1
#' @param impute_NA Logical value whether to impute NA values
#'
#' @importFrom rlang .data
#'
#' @keywords internal

impute_data <- function(aggregated_data, impute_to, impute_NA) {
  grouping_vars <- as.character(groups(aggregated_data))
  aggregated_data <- aggregated_data %>%
    group_by(.data$Metabolite) %>%
    mutate(imputed_Conc = zero_imputation(.data$Concentration, impute_to, impute_NA),
           .after = .data$Concentration) %>%
    group_by_at(grouping_vars)
  return(aggregated_data)
}


#' Transformation
#'
#' This function performs transformation of imputed concentration values
#' (imputed_Conc).
#'
#' @param vec MetAlyzer object
#' @param func A function for transformation
#'
#' @keywords internal

transform <- function(vec, func) {
  vec[vec > 0 & !is.na(vec)] <- func(vec[vec > 0 & !is.na(vec)])
  return(vec)
}


#' Transform aggregated data
#'
#' This function performs a transformation of imputed concentration values
#' (imputed_Conc) with a given function. NA values are skipped. The transformed
#' values are added to aggregated_data in an extra column transf_Conc.
#'
#' @param aggregated_data aggregated_data tibble data frame
#'
#' @importFrom rlang .data
#'
#' @keywords internal

transform_plotting_data <- function(aggregated_data) {
  aggregated_data <- mutate(aggregated_data,
                            log2_Conc = transform(.data$imputed_Conc, base::log2),
                            .after = .data$imputed_Conc)
  return(aggregated_data)
}
