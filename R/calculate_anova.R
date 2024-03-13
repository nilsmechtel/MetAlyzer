#' @title One-way ANOVA
#'
#' @description This method performs a one-way ANOVA on the grouped aggregated_data (the
#' categorical variable is removed from grouping first). The vector of the
#' categorical variable needs to have at least two levels after removing NAs
#' from the dependent variable vector. Otherwise a vector of NA is returned.
#' A Tukey post-hoc test is then used to determine group names, starting with
#' "A" followed by further letters. These group names are added to
#' aggregated_data in the column ANOVA_Group. Thereby, metabolites can be
#' identified which are significantly higher in one or more of the categorical
#' variable compared to all other for each metabolite.
#'
#' @param metalyzer_se A Metalyzer object
#' @param categorical A column defining the categorical variable
#' @param groups A vector of column names of aggregated_data to calculate the
#' ANOVA group wise. If the column does not exists in aggregated_data it is
#' automatically added from meta data. The default value is set to NULL, which
#' uses the existing grouping of aggregated_data.
#' @param impute_perc_of_min A numeric value below 1
#' @param impute_NA Logical value whether to impute NA values
#'
#' @return A data frame containing the log2 fold change for each metabolite
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
#' # reduced to only 'Acylcarnitines' (first metabolic class) for simplicity
#' drop_vec = unique(metalyzer_se@elementMetadata$metabolic_classes)[2:24]
#' metalyzer_se <- filterMetabolites(
#'   metalyzer_se,
#'   drop_metabolites = drop_vec
#' )
#' metalyzer_se <- filterMetaData(
#'   metalyzer_se,
#'   Tissue == "Drosophila"
#' )
#' metalyzer_se <- calculate_anova(
#'   metalyzer_se,
#'   categorical = "Extraction_Method",
#'   groups = c("Metabolite"),
#'   impute_perc_of_min = 0.2,
#'   impute_NA = TRUE
#' )
calculate_anova <- function(metalyzer_se, categorical, groups = NULL, impute_perc_of_min = 0.2, impute_NA = TRUE) {
  ## Check for agricolae installation
  # installed_packages <- utils::installed.packages()[, "Package"]
  # if (! "agricolae" %in% installed_packages) {
  #   cat("Installing package 'agricolae':\n")
  #   utils::install.packages("agricolae")
  #   cat("\n")
  # }
  
  ## Add grouping columns
  if (is.null(groups)) {
    groups <- as.character(groups(metalyzer_se@metadata$aggregated_data))
  } else {
    for (group in groups) {
      if (!group %in% colnames(metalyzer_se@metadata$aggregated_data)) {
        metalyzer_se <- expand_aggregated_data(metalyzer_se, meta_data_column = group)
      }
    }
  }

  ## Perform imputation and log2 transformation
  metalyzer_se <- data_imputation(metalyzer_se, impute_perc_of_min, impute_NA)
  metalyzer_se <- data_transformation(metalyzer_se)
  
  ## Add categorical column
  if (!categorical %in% colnames(metalyzer_se@metadata$aggregated_data)) {
    metalyzer_se <- expand_aggregated_data(metalyzer_se,
                                           meta_data_column = categorical)
  }
  
  ## Calculate ANOVA
  aggregated_data <- metalyzer_se@metadata$aggregated_data %>%
    group_by_at(unique(c(groups, categorical, "Metabolite"))) %>%
    mutate(ANOVA_n = sum(!is.na(.data$log2_Conc)))
  
  anova_data <- aggregated_data %>%
    dplyr::rename(Categorical = all_of(categorical)) %>%
    ungroup(.data$Categorical)
  cat(paste0("Info: Calculating ANOVA (groupwise: ",
             paste(groups(anova_data), collapse = " * "), ")...  "))
  anova_data <- mutate(
    anova_data,
    ANOVA_Group = calc_anova(.data$Categorical, .data$log2_Conc)
  )
  
  # c_vec <- filter(anova_data, Tissue == "Drosophila", Metabolite == "Cer(d18:1/23:0)")$Categorical
  # d_vec <- filter(anova_data, Tissue == "Drosophila", Metabolite == "Cer(d18:1/23:0)")$log2_Conc
  cat("finished!\n")

  if (any(aggregated_data$Concentration != anova_data$Concentration, na.rm = TRUE)) {
    stop("An unexpected error happened!")
  }
  aggregated_data$ANOVA_Group <- anova_data$ANOVA_Group

  metalyzer_se@metadata$aggregated_data <- aggregated_data
  return(metalyzer_se)
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
  ## if all concentration values equal to 0
  ## -> no imputation
  ## -> log2 transformation = NA
  ## -> no ANOVA is calculated (output: NA)
  
  # Check if all values in d_vec
  if (any(is.na(d_vec))) {
    # Return a placeholder value indicating ANOVA couldn't be performed
    return(NA)
  } else if (stats::var(d_vec) == 0) {
    # Check if all values are equal
    # Return a placeholder value indicating ANOVA couldn't be performed
    return(NA)
  } else if (length(unique(c_vec)) < 2) {
    # Check if c_vec has less than two levels
    # Return a placeholder value indicating ANOVA couldn't be performed
    return(NA)
  } else {
    tmp_df <- data.frame(Categorical = as.character(c_vec),
                         Dependent = as.numeric(d_vec))
    ## ANOVA
    anova <- stats::aov(Dependent ~ Categorical, data = tmp_df)
    ## Tukey post-hoc; each categorical variable gets assigned to a group
    anova_groups <- agricolae::HSD.test(anova, "Categorical", group = TRUE)$groups %>%
      select(-.data$Dependent) %>%
      tibble::rownames_to_column("Categorical") %>%
      mutate(Categorical = factor(.data$Categorical, levels = levels(c_vec))) %>%
      arrange(.data$Categorical) %>%
      tibble::deframe() %>%
      toupper()
    group_vec <- sapply(c_vec, function(m) anova_groups[m])
    return(group_vec)
  }
}
