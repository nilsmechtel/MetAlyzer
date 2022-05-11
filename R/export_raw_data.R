#' Export filtered raw data as csv
#'
#' This function exports the filtered raw data in the CSV format.
#'
#' @param object MetAlyzer object
#' @param sample_id A column from meta_data for sample identification
#' @param file_path file path
#'
#' @importFrom utils write.csv
#'
#' @keywords internal

export_raw_data <- function(object, sample_id, file_path) {
  meta_data <- metaData(object)
  raw_data <- rawData(object)
  df <- bind_cols(select(meta_data, all_of(sample_id)),
                  raw_data)
  cat("Number of samples:", nrow(meta_data), "\n")
  cat("Number of Metabolites:", ncol(raw_data), "\n")
  utils::write.csv(x = df,
                   file = file_path,
                   row.names = FALSE)
}
