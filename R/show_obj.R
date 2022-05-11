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
  meta_data <- metaData(object)
  quant_status <- quantStatus(object)
  metabolites <- metabolites(object)
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
  cat("\"MetAlyzer\" object:\n")
  cat("File name:", file, "\n")
  cat("Sheet:", sheet, "\n")
  cat("File path:", path, "\n")
  cat("Metabolites:", length(metabolites), "\n")
  cat("Classes:", length(unique(names(metabolites))), "\n")
  if (length(metabolites) > 0) {
    cat("Including metabolism indicators:",
        "Metabolism Indicators" %in% names(metabolites), "\n")
  }
  cat("Number of samples:", nrow(meta_data), "\n")
  if (ncol(meta_data) > 0) {
    cat(paste0("Columns meta data: \"", paste(colnames(meta_data),
                                              collapse = "\"; \""), "\"\n"))
  }
  cat("-------------------------------------\n")
}
