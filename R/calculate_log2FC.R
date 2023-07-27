#' Calculate log2 fold change
#'
#' This function calculates the log2 fold change of two groups from
#' plotting_data.
#' @param metalyzer A Metalyzer object
#' @param categorical A column specifying the two groups
#' @param installation_type A character, indicating the type of package to
#' @param impute_to A numeric value below 1
#' @param impute_NA Logical value whether to impute NA values
#' download and install. Options: ["binary", "source", "both"]
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
#' metalyzer_se <- MetAlyzer_dataset(file_path = extraction_data())
#' metalyzer_se <- renameMetaData(metalyzer_se, Method = `Sample Description`)
#' 
#' log2FC <- calculate_log2FC(metalyzer, Tissue, perc_of_min = 0.2, impute_NA = TRUE)

calculate_log2FC <- function(metalyzer_se, categorical, perc_of_min = 0.2, impute_NA = FALSE,
                             installation_type = "binary") {
  metalyzer_se <- impute_data(metalyzer_se, perc_of_min, impute_NA)
  metalyzer_se <- transform_plotting_data(metalyzer_se)
  aggregated_data <- metadata(metalyzer_se)$aggregated_data
  cat_str <- deparse(substitute(categorical))

  ## Check for qvalue and BiocManager installation
  installed_packages <- utils::installed.packages()[, "Package"]
  if (! "qvalue" %in% installed_packages) {
    if (! "BiocManager" %in% installed_packages) {
      cat("Installing package 'BiocManager':\n")
      utils::install.packages("BiocManager")
      cat("\n")
    }
    cat("Installing package 'qvalue':\n")
    BiocManager::install("qvalue", ask = FALSE, type = installation_type)
    cat("\n")
  }

  df <- aggregated_data %>%
    ungroup(all_of(cat_str)) %>%
    mutate(Value = log2_Conc,
        Group = !!sym(cat_str)) %>%
    select(Metabolite,
           Class,
           Group,
           Value)

  ## Check if already factor
  group_vec <- df$Group
  if (class(group_vec) != "factor") {
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
  cat("Calculating log2 fold change from ", group_levels[1], " to ",
      group_levels[2], " (column: ", cat_str, ").\n", sep = "")
  df <- filter(df, .data$Group %in% group_levels[1:2])
  df$Group <- droplevels(df$Group)

  ## Check for further grouping
  grouping_vars <- as.character(groups(df))

  if (!"Metabolite" %in% grouping_vars) {
    grouping_vars[length(grouping_vars)+1] <- "Metabolite"
  }
  cat("Calculating log2 fold change groupwise (",
      paste(grouping_vars, collapse = " * "), ").\n", sep = "")

  ## Calculate log2FC and p-val
  options(warn = -1)
  change_test_df <- df %>%
    group_modify(~ apply_linear_model(df = .x)) %>%
    ungroup(Metabolite) %>%
    mutate(qval = qvalue::qvalue(pval, pi0 = 1)$qvalues)
  options(warn = 0)
  return(change_test_df)
}


#' Calculate log2 fold change
#'
#' This function applies a linear model to calculate the log2 fold change and
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
    filter(!is.na(Value)) %>%
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