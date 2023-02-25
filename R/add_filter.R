#' Apply a filter on each group
#'
#' This function adds a new boolean column to aggregated_data based whether a
#' given criterion is met (TRUE) or not (FALSE).
#'
#' @param aggregated_data aggregated_data
#' @param ... Filter criterion, e.g. CV < 0.1
#'
#' @return An updated aggregated_data data frame
#'
#' @import dplyr
#' @export
#'
#' @examples
#' \dontrun{
#' print(1)
#' }

add_filter <- function(aggregated_data, ...) {
  orig_n_col <- ncol(aggregated_data)
  existing_filters <- grep("Filter[0-9]+", colnames(aggregated_data), value = TRUE)
  if (length(existing_filters) == 0) {
    new_filter_col <- "Filter1"
  } else {
    existing_i <- as.integer(stringr::str_extract(existing_filters, "[0-9]+"))
    new_filter_col <- paste0("Filter", max(existing_i, na.rm = TRUE) + 1)
  }
  cat("Adding new filter column:", new_filter_col, "\n")

  aggregated_data <- mutate(aggregated_data, ...)
  colnames(aggregated_data)[orig_n_col+1] <- new_filter_col
  return(aggregated_data)
}
