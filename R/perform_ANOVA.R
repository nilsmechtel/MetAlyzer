#' ANOVA
#'
#' This method performs a one-way ANOVA on the grouped plotting_data (the
#' categorical variable is removed from grouping first). For this, the column
#' valid_replicates must have at least one entry that is TRUE in each group.
#' Otherwise, a vector of NA is returned. A Tukey post-hoc test is then used to
#' determine group names, starting with "A" followed by further letters. These
#' group names are added to plotting_data in the column ANOVA_group. Thereby,
#' metabolites can be identified which are significantly higher in one or more
#' of the categorical variable compared to all other for each metabolite.
#'
#' @param object MetAlyzer object
#' @param categorical A column defining the categorical variable
#'
#' @importFrom rlang .data
#'
#' @keywords internal

perform_ANOVA <- function(object, categorical) {
  plotting_data <- object@plotting_data
  grouping_vars <- group_vars(plotting_data)
  plotting_data <- plotting_data %>%
    group_by_at(grouping_vars, add = FALSE) %>%
    ungroup(categorical) %>%
    mutate(ANOVA_group = calc_anova(get(categorical), .data$transf_Conc,
                                    .data$valid_replicates)) %>%
    group_by_at(grouping_vars)
  object@plotting_data <- plotting_data
  return(object)
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
#' @importFrom tibble rownames_to_column deframe
#' @importFrom stats aov
#' @importFrom agricolae HSD.test
#' @importFrom rlang .data
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
    groups <- HSD.test(anova, "Categorical", group = TRUE)$groups %>%
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
