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
