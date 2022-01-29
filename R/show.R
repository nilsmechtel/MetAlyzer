#' Show the 'MetAlyzer' object
#'
#' This function shows a summary of 'MetAlyzer' slot values.
#'
#' @param object MetAlyzer object
#'
#' @importFrom utils tail
#'
#' @keywords internal

show_obj <- function(object) {
  meta_data <- get_filtered_data(object, slot = "meta")
  quant_status <- get_filtered_data(object, slot = "quant")
  if (length(object@file_path) > 0) {
    s_fp <- strsplit(normalizePath(object@file_path), "/")[[1]]
    file <- tail(s_fp, 1)
    path <- paste(s_fp[1:length(s_fp)-1], collapse = "/")
  } else {
    file <- "<empty>"
    path <- "<empty>"
  }
  if (length(object@sheet) > 0) {
    sheet <- object@sheet
  } else {
    sheet <- "<empty>"
  }
  cat("-------------------------------------\n")
  cat("File name:", file, "\n")
  cat("Sheet:", sheet, "\n")
  cat("File path:", path, "\n")
  cat("Metabolites:", length(object@metabolites), "\n")
  cat("Classes:", length(unique(names(object@metabolites))), "\n")
  if (length(object@metabolites) > 0) {
    cat("Including metabolism indicators:", "Metabolism Indicators" %in% names(object@metabolites), "\n")
  }
  cat("Number of samples:", nrow(meta_data), "\n")
  if (ncol(meta_data) > 0) {
    cat(paste0("Columns meta data: \"", paste(colnames(meta_data), collapse = "\"; \""), "\"\n"))
  }
  cat("Ploting data created:", !nrow(object@plotting_data) == 0, "\n")
}
