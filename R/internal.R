#' Show the 'MetAlyzer' object
#'
#' This function shows a summary of 'MetAlyzer' slot values.
#'
#' @param object MetAlyzer object
#'
#' @importFrom utils tail
#'
#' @keywords internal

show_obj <- function(object) {
  meta_data <- metaData(object)
  quant_status <- quantStatus(object)
  metabolites <- metabolites(object)
  if (length(object@file_path) > 0) {
    s_fp <- strsplit(normalizePath(object@file_path), "/")[[1]]
    file <- tail(s_fp, 1)
    path <- paste(s_fp[1:length(s_fp)-1], collapse = "/")
  } else {
    file <- "<empty>"
    path <- "<empty>"
  }
  if (length(object@sheet) > 0) {
    sheet <- object@sheet
  } else {
    sheet <- "<empty>"
  }
  cat("-------------------------------------\n")
  cat("\"MetAlyzer\" object:\n")
  cat("File name:", file, "\n")
  cat("Sheet:", sheet, "\n")
  cat("File path:", path, "\n")
  cat("Metabolites:", length(metabolites), "\n")
  cat("Classes:", length(unique(names(metabolites))), "\n")
  if (length(metabolites) > 0) {
    cat("Including metabolism indicators:",
        "Metabolism Indicators" %in% names(metabolites), "\n")
  }
  cat("Number of samples:", nrow(meta_data), "\n")
  if (ncol(meta_data) > 0) {
    cat(paste0("Columns meta data: \"", paste(colnames(meta_data),
                                              collapse = "\"; \""), "\"\n"))
  }
  cat("-------------------------------------\n")
}


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


#' Print number of quantification status
#'
#' @param quant_status quant_status data frame
#' @param status quantification status
#' @param total total n
#'
#' @keywords internal

print_number <- function(quant_status, status, total) {
  number <- sum(quant_status == status, na.rm = TRUE)
  cat(paste0(status, ": ", number,
             " (", round(number/total*100, 2),"%)\n"))
}

#' Get all metabolites that have the quantification status at least once
#'
#' @param quant_status quant_status data frame
#' @param status quantification status
#'
#' @keywords internal

status_metabolites <- function(quant_status, status) {
  if (status == "NA") {
    n <- colSums(is.na(quant_status))
  } else {
    n <- colSums(quant_status == status)
  }
  metabolites <- colnames(quant_status)[n > 0]
  return(metabolites)
}

#' Summarize quantification status
#'
#' This function lists the number of each quantification status and its
#' percentage.
#'
#' @param object MetAlyzer object
#'
#' @keywords internal

sum_quant_data <- function(object) {
  quant_status <- quantStatus(object)
  nas <- sum(is.na(quant_status))
  total <- nrow(quant_status) * ncol(quant_status)
  status_vec <- levels(quant_status[,1])
  status_list <- list()
  cat("-------------------------------------\n")
  for (status in status_vec[which(status_vec %in% unlist(quant_status))]) {
    print_number(quant_status, status, total)
    status_list[[status]] <- status_metabolites(quant_status, status)
  }
  cat(paste0("NAs: ", nas, " (", round(nas/total*100, 2),"%)\n"))
  status_list[["NA"]] <- status_metabolites(quant_status, "NA")
  cat("-------------------------------------\n")
  return(status_list)
}


#' Filter meta data
#'
#' This function updates the "Filter" column in meta_data to filter out samples.
#'
#' @param object MetAlyzer object
#' @param column A length-one character vector specifying which column of
#' meta_data to use for filtering
#' @param keep A vector specifying which samples to keep
#' @param remove A vector specifying which samples to remove
#'
#' @keywords internal

filter_meta_data <- function(object, column, keep, remove) {
  old_filter <- object@meta_data$Filter
  if (!is.null(keep)) {
    new_filter <- object@meta_data[,column] %in% keep
  } else if (!is.null(remove)) {
    new_filter <- ! object@meta_data[,column] %in% remove
  }
  updated_filter <- sapply(1:nrow(object@meta_data), function(i) {
    old_filter[i] & new_filter[i]
  })
  object@meta_data$Filter <- updated_filter
  return(object)
}


#' Update meta data
#'
#' This function adds another column to filtered meta_data.
#'
#' @param object MetAlyzer object
#' @param name The new column name
#' @param new_colum A vector for the new column (length has to be same as the
#' number of filtered samples)
#'
#' @keywords internal

update_meta_data <- function(object, name, new_colum) {
  if (inherits(new_colum, "factor")) {
    levels <- levels(new_colum)
  } else {
    levels <- unique(new_colum)
  }
  meta_data <- object@meta_data
  meta_data[,name] <- factor(NA, levels = levels)
  meta_data[,name][meta_data$Filter == TRUE] <- new_colum
  object@meta_data <- meta_data
  return(object)
}


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
  metabolites <- object@metabolites[["filtered"]]
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
    cat("1 metabolite was removed!\n")
  } else if (diff > 1) {
    cat(paste(diff, "metabolites were removed!\n"))
  } else {
    cat("No metabolites were removed!\n")
  }
  object@metabolites[["filtered"]] <- metabolites
  return(object)
}


#' Reset metabolites
#'
#' This method resets the filtering of metabolites.
#'
#' @param object MetAlyzer object
#'
#' @keywords internal

reset_metabolites <- function(object) {
  original <- object@metabolites[["original"]]
  filtered <- object@metabolites[["filtered"]]
  diff <- length(original) - length(filtered)
  if (diff > 0) {
    cat(paste("Restoring", diff, "metabolite(s).\n"))
    object@metabolites[["filtered"]] <- original
  }
  return(object)
}


#' Get invalid metabolites
#'
#' This function extracts metabolites with a given percentage of valid replicates
#' from aggregated_data and filters metabolites.
#'
#' @param aggregated_data aggregated_data tibble data frame
#' @param filter_col Boolean column defining whether the filter criterion is met or not
#' @param replicate_t A numeric threshold
#'
#' @import dplyr
#' @importFrom rlang .data
#'
#' @keywords internal

get_invalid_metabolites <- function(aggregated_data, filter_col, replicate_t) {
  grouping_vars <- groups(aggregated_data)
  grouping_vars <- grouping_vars[-which(grouping_vars == "Metabolite")]

  cat(paste0("A metabolite counts as invalid, if ", round(replicate_t*100, 2),
             "% of replicates or less are considered as valid for it.\n"))
  if (length(grouping_vars) > 1) {
    cat("Replicates are split into each combination of grouping variables",
        paste0("(", paste(grouping_vars, collapse = " * "), ")\n"))
  }

  df <- aggregated_data %>%
    summarise(Filter_Col = first(!!sym(filter_col))) %>%
    group_by(.data$Metabolite) %>%
    summarise(n_valid = sum(.data$Filter_Col),
              n_tot = n(),
              perc = .data$n_valid/.data$n_tot) %>%
    filter(.data$perc <= replicate_t)
  invalid_metabolites <- as.character(df$Metabolite)
  return(invalid_metabolites)
}


#' Export filtered raw data as csv
#'
#' This function exports the filtered raw data in the CSV format.
#'
#' @param object MetAlyzer object
#' @param sample_id A column from meta_data for sample identification
#' @param file_path file path
#'
#' @importFrom utils write.csv
#'
#' @keywords internal

export_raw_data <- function(object, sample_id, file_path) {
  meta_data <- metaData(object)
  raw_data <- rawData(object)
  df <- bind_cols(select(meta_data, all_of(sample_id)),
                  raw_data)
  cat("Number of samples:", nrow(meta_data), "\n")
  cat("Number of Metabolites:", ncol(raw_data), "\n")
  utils::write.csv(x = df,
                   file = file_path,
                   row.names = FALSE)
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
