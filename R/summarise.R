#' Summarize quantification status
#'
#' This function lists the number of each quantification status and its
#' percentage.
#'
#' @param object MetAlyzer object
#'
#' @keywords internal

sum_quant_data <- function(object) {
  quant_status <- get_filtered_data(object, slot = "quant")
  nas <- sum(is.na(quant_status))
  total <- nrow(quant_status) * ncol(quant_status)
  print_number <- function(name) {
    number <- sum(quant_status == name, na.rm = TRUE)
    cat(paste0(name, ": ", number,
               " (", round(number/total*100, 2),"%)\n"))
  }
  cat("-------------------------------------\n")
  status_vec <- levels(quant_status[,1])
  for (status in status_vec[which(status_vec %in% unlist(quant_status))]) {
    print_number(status)
  }
  cat(paste0("NAs: ", nas, " (", round(nas/total*100, 2),"%)\n"))
}
