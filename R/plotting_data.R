#' Reshape data
#'
#' This function reshapes raw_data, quant_status and meta_data and combines them
#' in a tibble data frame for plotting with ggplot2.
#' @param object MetAlyzer object
#' @param ... A selection of columns from meta_data to add to reshaped data frame
#'
#' @import dplyr
#' @import tidyr
#'
#' @keywords internal

plotting_data <- function(object, ...) {
  meta_data <- get_filtered_data(object, slot = "meta", verbose = FALSE)
  raw_data <- get_filtered_data(object, slot = "data", verbose = FALSE)
  quant_status <- get_filtered_data(object, slot = "quant", verbose = FALSE)
  if (all(nrow(raw_data)>0, nrow(quant_status)>0)) {
    extra_columns <- select(meta_data, ...)
    raw_data <- bind_cols(raw_data, extra_columns)
    gathered_data <- gather(raw_data, key = Metabolite, value = Concentration, -c(...))
    gathered_status <- gather(quant_status, key = Metabolite, value = Status)
    col_names <- colnames(gathered_data)
    group_cols <- col_names[1:which(col_names == "Metabolite")]
    plotting_data <- gathered_data %>%
      group_by_at(group_cols) %>%
      mutate(Class = sapply(Metabolite, function(x) {
        names(object@metabolites[object@metabolites == x])
      }),
      Replicates = 1:n(),
      .after = Metabolite)
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
      } else {
        levels <- unique(col)
      }
      plotting_data[, col_name] <- factor(unlist(plotting_data[,col_name]),
                                          levels = levels)
    }
    return(plotting_data)
  }
}


#' Threshold CV
#'
#' This function assigns a CV value according to a vector of thresholds.
#' @param x A CV value
#' @param ts A numeric vector of thresholds
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
#'
#' @keywords internal

calc_CV <- function(plotting_data, ts) {
  plotting_data <- plotting_data %>%
    mutate(Mean = mean(Concentration, na.rm = TRUE),
           SD = sd(Concentration, na.rm = TRUE),
           CV = SD / Mean,
           CV_thresh = sapply(CV, function(x) set_threshold(x, ts = ts)),
           .after = Concentration)
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
#'
#' @keywords internal

valid_measurement <- function(plotting_data, valid_vec, t) {
  filter_status <- function(vec) {
    sum(vec %in% valid_vec) / length(vec) > t
  }
  plotting_data <- mutate(plotting_data,
                          Valid = filter_status(Status),
                          .after = Status)
  return(plotting_data)
}


#' Zero imputation
#'
#' This function performs zero imputation with the minimal positive value times i
#' @param vec A numeric vector containing the concentration values
#' @param i A numeric value below 1)
#' @param imputeNA Logical value whether to impute NA values. Default = FALSE
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


#' Logarithmic transformation
#'
#' This function performs logarithmic transformation of imputed concentration
#' values (imp_Conc)
#' @param vec MetAlyzer object
#' @param func A logarithmic function
#'
#' @keywords internal

log_transform <- function(vec, func) {
  vec[vec > 0 & !is.na(vec)] <- func(vec[vec > 0 & !is.na(vec)])
  return(vec)
}


#' Perform an ANOVA
#'
#' This function filters based on the filter vector valid_vec, performs a
#' one-way ANOVA and adds the column Group to plotting_data with the results of
#' a Tukey post-hoc test
#' @param c_vec A character vector containing the categorical variables
#' @param d_vec A numeric vector containing the dependent variables
#' @param valid_vec A logical vector for filtering
#'
#' @import dplyr
#' @import tibble
#'
#' @keywords internal

calc_anova <- function(c_vec, d_vec, valid_vec) {
  # if all concentration values equal 0 (no imputation; log -> NA) or no method
  # achieves valid concentrations no ANOVA is calculated (output: NA)
  if (all(is.na(d_vec)) | sum(valid_vec) == 0) {
    group_vec <- as.character(rep(NA, length(c_vec)))
  } else {
    tmp_df <- data.frame(Categorical = as.character(c_vec),
                         Dependent = as.numeric(d_vec))
    # # ANOVA
    anova <- aov(Dependent ~ Categorical, data = tmp_df)
    # Tukey post-hoc; each categorical variable gets assigned to a group
    groups <- agricolae::HSD.test(anova, "Categorical", group = TRUE)$groups %>%
      select(-Dependent) %>%
      rownames_to_column("Categorical") %>%
      mutate(Categorical = factor(Categorical, levels = levels(c_vec))) %>%
      arrange(Categorical) %>%
      deframe() %>%
      toupper()
    group_vec <- sapply(c_vec, function(m) groups[m])
  }
  return(group_vec)
}
