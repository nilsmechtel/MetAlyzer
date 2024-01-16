#' @title Add column to aggregated_data
#'
#' @description This function adds a column to the aggregated_data tibble by
#' using information from the meta data. 
#'
#' @param metalyzer_se A Metalyzer object
#' @param meta_data_column A column from meta data
#' @return An updated Metalyzer object
#' 
#' @import dplyr
#' @importFrom rlang .data
#' 
#' @keywords internal
expand_aggregated_data <- function(metalyzer_se, meta_data_column) {
  aggregated_data <- metalyzer_se@metadata$aggregated_data
  meta_data <- colData(metalyzer_se)
  if (!meta_data_column %in% colnames(meta_data)) {
    cat("Warning: Could not find column", meta_data_column, "in meta data!\n")
  } else if (meta_data_column %in% colnames(aggregated_data)) {
    cat("Info: Column", meta_data_column, "already exists in aggregated_data.\n")
  } else {
    if (!is.character(meta_data_column)) {
      meta_data_column <- deparse(substitute(meta_data_column))
    }
    mapping_vec <- unlist(meta_data[meta_data_column])
    names(mapping_vec) <- rownames(meta_data[meta_data_column])
    aggregated_data <- dplyr::mutate(aggregated_data, 
                                     !!meta_data_column := factor(sapply(.data$ID, function(id) {
                                       mapping_vec[id]
                                     }),
                                     levels = unique(mapping_vec)), .after = .data$ID)
    
    metalyzer_se@metadata$aggregated_data <- aggregated_data
  }
  return(metalyzer_se)
}
