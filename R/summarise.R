#' Show a MetAlyzer object
#'
#' This function shows a summary of MetAlyzer slot values
#' @param object MetAlyzer object
#'
#' @return
#' @export
#'
#' @examples

show_obj <- function(object) {
  if (length(object@file_path) > 0) {
    s_fp <- strsplit(object@file_path, "/")[[1]]
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
  cat("Number of samples:", nrow(object@meta_data), "\n")
  if (ncol(object@meta_data) > 0) {
    cat("Columns meta data:", paste(colnames(object@meta_data), collapse = "; "), "\n")
  }
  cat("Quantification status completed:", !nrow(object@quant_status) == 0, "\n")
  cat("-------------------------------------\n")
}


#' Summarize quantification status
#'
#' This function lists the proportion of LODs, LOQs, Valids, calibration range passes and NAs
#' @param object MetAlyzer object
#'
#' @return
#' @export
#'
#' @examples

sum_quant_data <- function(object) {
  sum_LOD <- sum(object@quant_status == "LOD", na.rm = TRUE)
  sum_LOQ <- sum(object@quant_status == "LOQ", na.rm = TRUE)
  sum_valid <- sum(object@quant_status == "Valid", na.rm = TRUE)
  sum_val_range <- sum(object@quant_status == "Out of calibration range", na.rm = TRUE)
  nas <- sum(is.na(object@quant_status))
  total <- nrow(object@quant_status)*ncol(object@quant_status)
  cat("-------------------------------------\n")
  cat(paste0("Number of LODs: ", sum_LOD,
               " (", round(sum_LOD/total*100, 2), "%)\n"))
  cat(paste0("Number of LOQs: ", sum_LOQ,
               " (", round(sum_LOQ/total*100, 2),"%)\n"))
  cat(paste0("Number of Valids: ", sum_valid,
               " (", round(sum_valid/total*100, 2),"%)\n"))
  cat(paste0("Number of calibration range passes: ", sum_val_range,
               " (", round(sum_val_range/total*100, 2),"%)\n"))
  cat(paste0("Number of NAs: ", nas,
               " (", round(nas/total*100, 2),"%)\n"))
  cat("-------------------------------------\n")
}
