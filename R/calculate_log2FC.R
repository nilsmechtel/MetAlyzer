#' Calculate log2 fold change
#'
#' This function calculates the log2 fold change of two groups from
#' plotting_data.
#' @param aggregated_data aggregated_data
#' @param categorical A column specifying the two groups
#' @param installation_type A character, indicating the type of package to
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
#' \dontrun{
#' print(1)
#' }

calculate_log2FC <- function(aggregated_data, categorical,
                             installation_type = "binary") {
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
    rename(Value = .data$log2_Conc,
           Group = !!sym(cat_str)) %>%
    select(.data$Metabolite,
           .data$Class,
           .data$Group,
           .data$Value)

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
    ungroup(.data$Metabolite) %>%
    mutate(qval = qvalue::qvalue(.data$pval, pi0 = 1)$qvalues)
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
