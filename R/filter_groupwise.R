#' Apply a filter on each group
#'
#' This function takes a boolean column from aggregated_data and filters based
#' on the percentage of samples within each group that meet the filter criterion
#' (=TRUE).
#'
#' @param aggregated_data aggregated_data
#' @param ... Columns of aggregated_data to group by
#' @param filter_col Column with logical values that defines whether a filter is
#' met or not
#' @param lower_bound Numerical threshold between 0 and 1 defining the lower
#' bound (t < x)
#' @param upper_bound Numerical threshold between 0 and 1 defining the upper
#' bound (x < t)
#'
#' @return A filtered aggregated_data data frame
#'
#' @import dplyr
#' @importFrom rlang .data
#' @export

filter_groupwise <- function(aggregated_data, ...,
                             filter_col=.data$Valid_Replicates,
                             lower_bound=0,
                             upper_bound=FALSE) {
  col_str <- deparse(substitute(filter_col))
  grouping_vars <- as.character(groups(aggregated_data))

  aggregated_data <- aggregated_data %>%
    group_by(...) %>%
    mutate(filter_col = sum(!!sym(col_str))/n()) %>%
    group_by_at(grouping_vars)

  if (!isFALSE(lower_bound)) {
    aggregated_data <- filter(aggregated_data, .data$filter_col > lower_bound)
  }
  if (!isFALSE(upper_bound)) {
    aggregated_data <- filter(aggregated_data, .data$filter_col < upper_bound)
  }

  aggregated_data <- select(aggregated_data, -.data$filter_col)
  return(aggregated_data)
}
