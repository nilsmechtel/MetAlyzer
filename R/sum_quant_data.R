#' Summarize quantification status
#'
#' This function lists the number of each quantification status and its
#' percentage.
#'
#' @param object MetAlyzer object
#'
#' @keywords internal

sum_quant_data <- function(object) {
  quant_status <- quantStatus(object)
  nas <- sum(is.na(quant_status))
  total <- nrow(quant_status) * ncol(quant_status)
  status_vec <- levels(quant_status[,1])
  status_list <- list()
  cat("-------------------------------------\n")
  for (status in status_vec[which(status_vec %in% unlist(quant_status))]) {
    print_number(quant_status, status, total)
    status_list[[status]] <- status_metabolites(quant_status, status)
  }
  cat(paste0("NAs: ", nas, " (", round(nas/total*100, 2),"%)\n"))
  status_list[["NA"]] <- status_metabolites(quant_status, "NA")
  cat("-------------------------------------\n")
  return(status_list)
}


#' Print number of quantification status
#'
#' @param quant_status quant_status data frame
#' @param status quantification status
#' @param total total n
#'
#' @keywords internal

print_number <- function(quant_status, status, total) {
  number <- sum(quant_status == status, na.rm = TRUE)
  cat(paste0(status, ": ", number,
             " (", round(number/total*100, 2),"%)\n"))
}


#' Get all metabolites that have the quantification status at least once
#'
#' @param quant_status quant_status data frame
#' @param status quantification status
#'
#' @keywords internal

status_metabolites <- function(quant_status, status) {
  if (status == "NA") {
    n <- colSums(is.na(quant_status))
  } else {
    n <- colSums(quant_status == status)
  }
  metabolites <- colnames(quant_status)[n > 0]
  return(metabolites)
}


