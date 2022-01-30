#' Create plotting data
#'
#' This function reshapes raw_data, quant_status and meta_data and combines them
#' together with basic statistics in a tibble data frame for plotting with
#' 'ggplot2'. plotting_data is grouped by metabolites as well as the selection of
#' additional variables. Statistics are then calculated per group.
#' @param object MetAlyzer object
#' @param ... A selection of columns from meta_data to add to reshaped data frame
#' @param ungrouped A column from meta_data to add to reshaped data frame that
#' will not be used as grouping variables
#' @param ts A numeric vector of thresholds for CV categorization
#' @param valid_vec A character vector containing each quantification status that
#' is considered to be a valid measurement
#' @param t A numeric threshold to determine valid measurements
#'
#' @keywords internal

create_plotting_data <- function(object, ..., ungrouped, ts, valid_vec, t) {
  plotting_data <- plotting_data(object, ..., ungrouped = ungrouped)
  plotting_data <- calc_CV(plotting_data, ts = ts)
  plotting_data <- valid_measurement(plotting_data, valid_vec, t)
  object@plotting_data <- plotting_data
  return(object)
}


#' Reshape data
#'
#' This function reshapes raw_data, quant_status and meta_data and combines them
#' in a tibble data frame for plotting with 'ggplot2'.
#' @param object MetAlyzer object
#' @param ... A selection of columns from meta_data to add to reshaped data frame
#' @param ungrouped A column from meta_data to add to reshaped data frame that
#' will not be used as grouping variables
#'
#' @import dplyr
#' @importFrom tidyr gather
#' @importFrom rlang .data
#'
#' @keywords internal

plotting_data <- function(object, ..., ungrouped) {
  meta_data <- get_filtered_data(object, slot = "meta")
  raw_data <- get_filtered_data(object, slot = "data")
  quant_status <- get_filtered_data(object, slot = "quant")
  if (nchar(ungrouped)==0 | ungrouped == "NULL") {ungrouped <- NULL}
  if (nrow(meta_data) > 0) {
    extra_columns <- select(meta_data, c(... , ungrouped))
    raw_data <- bind_cols(raw_data, extra_columns)
    gathered_data <- gather(raw_data, key = "Metabolite",
                            value = "Concentration", -colnames(extra_columns))
    gathered_status <- gather(quant_status, key = "Metabolite",
                              value = "Status")
    col_names <- colnames(select(gathered_data, -ungrouped))
    group_cols <- col_names[1:which(col_names == "Metabolite")]
    plotting_data <- gathered_data %>%
      group_by_at(group_cols) %>%
      mutate(Class = sapply(.data$Metabolite, function(x) {
        names(object@metabolites[object@metabolites == x])
      }),
      .after = .data$Metabolite)
    plotting_data$Metabolite <- factor(plotting_data$Metabolite,
                                       levels = unique(object@metabolites))
    plotting_data$Class <- factor(plotting_data$Class,
                                  levels = unique(names(object@metabolites)))
    plotting_data$Status <- factor(gathered_status$Status,
                                   levels = levels(quant_status[,1]))
    for (col_name in colnames(extra_columns)) {
      col <- extra_columns[, col_name]
      if (class(col) == "factor") {
        levels <- levels(col)
      } else if (!any(grepl("\\D", col))) {
        levels <- unique(as.character(sort(as.numeric(col))))
      } else {
        levels <- unique(col)
      }
      plotting_data[, col_name] <- factor(unlist(plotting_data[,col_name]),
                                          levels = levels)
    }
    plotting_data <- arrange_at(plotting_data, c(group_cols, ungrouped))

    return(plotting_data)
  } else {
    return(data.frame())
  }
}


#' Threshold CV
#'
#' This function assigns a CV value according to a vector of thresholds.
#' @param x A CV value
#' @param ts A numeric vector of thresholds
#'
#' @importFrom utils tail
#'
#' @keywords internal

set_threshold <- function(x, ts) {
  levels <- c(paste0("max", ts*100), paste0("more", tail(ts, 1)*100))
  ts <- sort(ts, decreasing = TRUE)
  if (is.na(x)) {
    v <- NA
  } else {
    for (t in ts) {
      if (x <= t) {
        v <- paste0("max", t*100)
      }
    }
    if (x > ts[1]) {
      v <- paste0("more", ts[1]*100)
    }
  }
  v <- factor(v, levels = levels)
  return(v)
}


#' Add CV
#'
#' This function adds the mean, standard deviation (SD) and the
#' coefficient of variation (CV) to plotting_data
#' @param plotting_data plotting_data tibble data frame
#' @param ts A numeric vector of thresholds
#'
#' @import dplyr
#' @importFrom stats sd
#' @importFrom rlang .data
#'
#' @keywords internal

calc_CV <- function(plotting_data, ts) {
  plotting_data <- plotting_data %>%
    mutate(Mean = mean(.data$Concentration, na.rm = TRUE),
           SD = sd(.data$Concentration, na.rm = TRUE),
           CV = .data$SD / .data$Mean,
           CV_thresh = sapply(.data$CV, function(x) set_threshold(x, ts = ts)),
           .after = .data$Concentration)
  return(plotting_data)
}


#' Add valid filter
#'
#' This function adds a filter column based on the quantification status. The
#' filter is True if the percentage of measurements with a quantification status
#' part of valid_vec is greater than the threshold t.
#' @param plotting_data plotting_data tibble data frame
#' @param valid_vec A character vector containing each quantification status that
#' is considered to be a valid measurement
#' @param t A numeric threshold
#'
#' @import dplyr
#' @importFrom rlang .data
#'
#' @keywords internal

valid_measurement <- function(plotting_data, valid_vec, t) {
  filter_status <- function(vec) {
    sum(vec %in% valid_vec) / length(vec) > t
  }
  plotting_data <- mutate(plotting_data,
                          valid_replicates = filter_status(.data$Status),
                          .after = .data$Status)
  return(plotting_data)
}
