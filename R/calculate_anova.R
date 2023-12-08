#' One-way ANOVA
#'
#' This method performs a one-way ANOVA on the grouped aggregated_data (the
#' categorical variable is removed from grouping first). The vector of the
#' categorical variable needs to have at least two levels after removing NAs
#' from the dependent variable vector. Otherwise a vector of NA is returned.
#' A Tukey post-hoc test is then used to determine group names, starting with
#' "A" followed by further letters. These group names are added to
#' aggregated_data in the column ANOVA_Group. Thereby, metabolites can be
#' identified which are significantly higher in one or more of the categorical
#' variable compared to all other for each metabolite.
#'
#' @param aggregated_data aggregated_data
#' @param categorical A column defining the categorical variable
#'
#' @return A data frame containing the log2 fold change for each metabolite
#'
#' @import dplyr
#' @importFrom utils install.packages
#' @importFrom utils installed.packages
#' @importFrom rlang .data
#' @export
#' 
#' @examples
#' metalyzer_se <- MetAlyzer_dataset(file_path = example_extraction_data())
#! examples
calculate_anova <- function(aggregated_data, categorical, impute_perc_of_min = 0.2, impute_NA = FALSE) {
  cat_str <- deparse(substitute(categorical))

  ## Check for agricolae installation
  installed_packages <- utils::installed.packages()[, "Package"]
  if (! "agricolae" %in% installed_packages) {
    cat("Installing package 'agricolae':\n")
    utils::install.packages("agricolae")
    cat("\n")
  }

  # metalyzer_se <- data_imputation(metalyzer_se, impute_perc_of_min, impute_NA)
  # metalyzer_se <- data_transformation(metalyzer_se)

  anova_data <- aggregated_data %>%
    rename(Categorical = all_of(cat_str)) %>%
    ungroup(.data$Categorical) %>%
    mutate(ANOVA_Group = calc_anova(.data$Categorical,
                                    .data$log2_Conc))

  if (any(aggregated_data$Concentration != anova_data$Concentration, na.rm = TRUE)) {
    stop("An unexpected error happened!")
  }
  aggregated_data$ANOVA_Group <- anova_data$ANOVA_Group

  return(aggregated_data)
}


#' Perform an ANOVA
#'
#' This function filters based on the filter vector valid_vec, performs a
#' one-way ANOVA and adds the column Group to aggregated_data with the results of
#' a Tukey post-hoc test
#' @param c_vec A character vector containing the categorical variables
#' @param d_vec A numeric vector containing the dependent variables
#'
#' @import dplyr
#' @importFrom tibble rownames_to_column deframe
#' @importFrom stats aov
#' @importFrom rlang .data
#'
#' @keywords internal

calc_anova <- function(c_vec, d_vec) {
  ## if all concentration values equal to 0 (no imputation; log -> NA) no ANOVA
  ## is calculated (output: NA)

  ## filter out NA
  c_vec <- c_vec[!is.na(d_vec)]
  d_vec <- d_vec[!is.na(d_vec)]

  if (length(unique(c_vec)) < 2) {
    ## less than two levels
    group_vec <- as.character(rep(NA, length(c_vec)))
  } else {
    tmp_df <- data.frame(Categorical = as.character(c_vec),
                         Dependent = as.numeric(d_vec))
    ## ANOVA
    anova <- aov(Dependent ~ Categorical, data = tmp_df)
    ## Tukey post-hoc; each categorical variable gets assigned to a group
    groups <- agricolae::HSD.test(anova, "Categorical", group = TRUE)$groups %>%
      select(-.data$Dependent) %>%
      rownames_to_column("Categorical") %>%
      mutate(Categorical = factor(.data$Categorical, levels = levels(c_vec))) %>%
      arrange(.data$Categorical) %>%
      deframe() %>%
      toupper()
    group_vec <- sapply(c_vec, function(m) groups[m])
  }
  return(group_vec)
}
