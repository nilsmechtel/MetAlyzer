#' Reshape data
#'
#' This function reshapes raw_data, quant_status and meta_data and combines them
#' in a tibble data frame for plotting with ggplot2.
#' @param object MetAlyzer object
#' @param ... A selection of columns from meta_data to add to reshaped data frame
#'
#' @return
#'
#' @import dplyr
#' @import tidyr
#' @export

plotting_data <- function(object, ...) {
  meta_data <- get_filtered_data(object, slot = "meta", verbose = FALSE)
  raw_data <- get_filtered_data(object, slot = "data", verbose = FALSE)
  quant_status <- get_filtered_data(object, slot = "quant", verbose = FALSE)
  if (all(nrow(raw_data)>0, nrow(quant_status)>0)) {
    extra_columns <- select(meta_data, ...)
    raw_data <- bind_cols(raw_data, extra_columns)
    gathered_data <- gather(raw_data, key = Metabolite, value = Concentration, -c(...))
    gathered_status <- gather(quant_status, key = Metabolite, value = Status)
    plotting_data <- mutate(gathered_data,
                            Class = sapply(Metabolite, function(x) {
                              names(object@metabolites[object@metabolites == x])
                            }),
                            .after = Metabolite)
    plotting_data$Metabolite <- factor(plotting_data$Metabolite, levels = unique(object@metabolites))
    plotting_data$Class <- factor(plotting_data$Class, levels = unique(names(object@metabolites)))
    plotting_data$Status <- factor(gathered_status$Status, levels = unique(gathered_status$Status))
    return(plotting_data)
  }
}


#' Threshold CV
#'
#' This function assigns a CV value according to a vector of thresholds.
#' @param x A CV value
#' @param ts A numeric vector of thresholds
#'
#' @return
#'
#' @export

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
#' @return
#'
#' @import dplyr
#' @export

calc_CV <- function(plotting_data, ts) {
  col_names <- colnames(plotting_data)
  group_cols <- col_names[1:which(col_names == "Metabolite")]
  plotting_data <- plotting_data %>%
    group_by_at(group_cols) %>%
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
#' @return
#'
#' @import dplyr
#' @export

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
#'
#' @return
#'
#' @export

zero_imputation <- function(vec, i) {
  non_zero <- vec[vec > 0]
  imp_v <- ifelse(length(non_zero) > 0, min(non_zero) * i, NA)
  vec[vec == 0] <- imp_v
  return(vec)
}


#' Logarithmic transformation
#'
#' This function performs logarithmic transformation of imputed concentration
#' values (imp_Conc)
#' @param vec MetAlyzer object
#' @param func A logarithmic function
#'
#' @return
#'
#' @export

log_transform <- function(vec, func) {
  vec[vec > 0 & !is.na(vec)] <- func(vec[vec > 0 & !is.na(vec)])
  return(vec)
}


#' Perform an ANOVA
#'
#' This function filters based on the filter vector valid_vec and then performs
#' a one-way ANOVA
#' @param c_vec A character vector containing the categorical variables
#' @param d_vec A numeric vector containing the dependent variables
#' @param valid_vec A logical vector for filtering
#'
#' @return
#'
#' @import dplyr
#' @import tibble
#' @export

calc_anova <- function(c_vec, d_vec, valid_vec) {
  # if all concentration values equal 0 (no imputation; log -> NA) or no method
  # achieves valid concentrations no ANOVA is calculated (output: NA)
  if (all(is.na(d_vec)) | sum(valid_vec) == 0) {
    group_vec <- as.character(rep(NA, length(c_vec)))
  } else {
    tmp_df <- data.frame(Categorical = as.character(c_vec),
                         Dependent = as.numeric(d_vec))
    # ANOVA
    anova <- aov(Dependent ~ Categorical, data = tmp_df)
    # Tukey post-hoc; each categorical variable gets assigned to a group
    groups <- agricolae::HSD.test(anova, "Categorical", group=TRUE)$groups %>%
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
