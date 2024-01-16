#' @title Calculate log2 fold change
#'
#' @description This function calculates the log2 fold change of two groups from
#' plotting_data.
#' @param metalyzer_se A Metalyzer object
#' @param categorical A column specifying the two groups
#' @param impute_perc_of_min A numeric value below 1
#' @param impute_NA Logical value whether to impute NA values
#'
#' @return A data frame containing the log2 fold change for each metabolite
#'
#' @import dplyr
#' @import SummarizedExperiment
#' @importFrom rlang .data
#' @importFrom data.table :=
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer_dataset(file_path = example_mutation_data_xl())
#' metalyzer_se <- filterMetabolites(
#'   metalyzer_se,
#'   drop_metabolites = "Metabolism Indicators"
#' )
#' metalyzer_se <- renameMetaData(
#'   metalyzer_se,
#'   Mutant_Control = "Sample Description"
#' )
#' 
#' metalyzer_se <- calculate_log2FC(
#'   metalyzer_se,
#'   categorical = "Mutant_Control",
#'   impute_perc_of_min = 0.2,
#'   impute_NA = FALSE
#' )
calculate_log2FC <- function(metalyzer_se,
                             categorical,
                             impute_perc_of_min = 0.2,
                             impute_NA = FALSE) {
  ## Check for qvalue and BiocManager installation
  # installed_packages <- utils::installed.packages()[, "Package"]
  # if (! "qvalue" %in% installed_packages) {
  #   if (! "BiocManager" %in% installed_packages) {
  #     cat("Info: Installing package 'BiocManager':\n")
  #     utils::install.packages("BiocManager")
  #     cat("\n")
  #   }
  #   cat("Info: Installing package 'qvalue':\n")
  #   BiocManager::install("qvalue", ask = FALSE, type = "binary")
  #   cat("\n")
  # }
  
  
  ## Create a new dataframe to calculate the log2FC
  if (!categorical %in% colnames(metalyzer_se@metadata$aggregated_data)) {
    metalyzer_se <- expand_aggregated_data(metalyzer_se,
                                           meta_data_column = categorical)
  }
  
  metalyzer_se <- data_imputation(metalyzer_se, impute_perc_of_min, impute_NA)
  metalyzer_se <- data_transformation(metalyzer_se)

  df <- metalyzer_se@metadata$aggregated_data %>%
    ungroup(all_of(categorical)) %>%
    mutate(Value = .data$log2_Conc,
        Group = !!sym(categorical)) %>%
    select(.data$Metabolite,
           .data$Class,
           .data$Group,
           .data$Value)

  ## Check if already factor
  group_vec <- df$Group
  if (!methods::is(group_vec, "factor")) {
    group_vec <- factor(group_vec)
    df$Group <- group_vec
    cat("Warning: No order was given for categorical!\n")
  } else {
    group_vec <- droplevels(group_vec)
  }
  group_levels <- levels(group_vec)
  if (length(group_levels) > 2) {
    cat("Warning: More than two levels were given! Dropping",
        paste(group_levels[3:length(group_levels)], collapse = ", "), "\n")
  }
  cat("Info: Calculating log2 fold change from ", group_levels[1], " to ",
      group_levels[2], " (column: ", categorical, ").\n", sep = "")
  df <- filter(df, .data$Group %in% group_levels[1:2])
  df$Group <- droplevels(df$Group)

  ## Check for further grouping
  grouping_vars <- as.character(groups(df))

  if (!"Metabolite" %in% grouping_vars) {
    grouping_vars[length(grouping_vars)+1] <- "Metabolite"
  }
  
  ## Calculate log2FC and p-val
  cat("Info: Calculating log2 fold change groupwise (",
      paste(grouping_vars, collapse = " * "),
      ") using a linear model...  ", sep = "")
  options(warn = -1)
  change_df <- df %>%
    group_modify(~ apply_linear_model(df = .x)) %>%
    ungroup(.data$Metabolite) %>%
    mutate(qval = qvalue::qvalue(.data$pval, pi0 = 1)$qvalues)
  options(warn = 0)
  cat("finished!\n")
  
  metalyzer_se@metadata$log2FC <- change_df
  return(metalyzer_se)
}

#' @title Calculate log2 fold change
#'
#' @description This function applies a linear model to calculate the log2 fold change and
#' its significance
#' @param df A subset data frame
#' @param ...
#'
#' @import dplyr
#' @importFrom stats lm
#' @importFrom rlang .data
#'
#' @keywords internal
apply_linear_model <- function(df, ...) {
  class <- df$Class[1]
  df <- df %>%
    filter(!is.na(.data$Value)) %>%
    droplevels()
  if (length(levels(df$Group)) != 2) {
    l2fc <- NA
    pval <- NA
  } else {
    fit1 <- stats::lm(Value ~ Group, data = df)
    l2fc <- fit1$coefficients[2]
    fit_dim <- dim(summary(fit1)$coefficients)
    pval <- summary(fit1)$coefficients[fit_dim[1],fit_dim[2]]
  }
  output_df <- data.frame(Class = class,
                          log2FC = l2fc,
                          pval = pval,
                          row.names = NULL)
  return(output_df)
}