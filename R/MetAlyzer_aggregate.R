#' @title Reshape data
#'
#' @description This function reshapes selected meta_data variables,
#' conc_values, quant_status and metatabolites and combines them into a
#' tibble data frame for filtering with dplyr and plotting with 'ggplot2'.
#'
#' @param metabolites metabolites MetAlyzer object
#' @param meta_columns A selection of columns from meta_data to add to
#' aggregated data frame
#' @param conc_values conc_values of a MetAlyzer object
#' @param quant_status quant_status of a MetAlyzer object
#' @import dplyr
#'
#' @keywords internal

reshape_data <- function(metabolites,
                        classes,
                        meta_columns,
                        conc_values,
                        quant_status) {
    group_cols <- c(colnames(meta_columns), "Metabolite")
    meta_columns <- dplyr::mutate(
      meta_columns,
      ID = rownames(meta_columns),
      .before = 1
    )
    comb_data <- dplyr::bind_cols(meta_columns, conc_values)
    gathered_data <- tidyr::gather(comb_data, key = "Metabolite",
                            value = "Concentration", -colnames(meta_columns))
    gathered_status <- tidyr::gather(quant_status, key = "Metabolite",
                              value = "Status")

    aggregated_data <- gathered_data %>%
      dplyr::group_by_at(group_cols) %>%
      dplyr::mutate(Class = sapply(Metabolite, function(x) {
        classes[metabolites == x]
      }),
      .after = Metabolite)

    aggregated_data$Metabolite <- factor(aggregated_data$Metabolite,
                                          levels = unique(metabolites))
    aggregated_data$Class <- factor(aggregated_data$Class,
                                    levels = unique(classes))
    aggregated_data$Status <- factor(gathered_status$Status,
                                      levels = levels(quant_status[,1]))
    for (col_name in colnames(meta_columns)) {
      col <- meta_columns[, col_name]
      if (inherits(col, "factor")) {
        levels <- levels(col)
      } else if (!any(grepl("\\D", col))) {
        levels <- unique(as.character(sort(as.numeric(col))))
      } else {
        levels <- unique(col)
      }
      aggregated_data[, col_name] <- factor(unlist(aggregated_data[,col_name]),
                                            levels = levels)
    }
    aggregated_data <- arrange_at(aggregated_data, group_cols)
  
  return(droplevels(aggregated_data))
}


#' @title Add valid filter
#'
#' @description This function adds a filter column to aggregated_data based
#' on the quantification status. The filter is True if the percentage of
#' valid measurements is at least a given threshold.
#'
#' @param aggregated_data aggregated_data tibble data frame
#' @param status_list A list with three entries.
#' "Valid" holds the vector with a quantification status that is valid.
#' "Valid threshold" holds the minimal percentage threshold.
#' "KO" holds a character vector with quantification status that directly
#' invalidates a set of repetitions once it is encountered.
#' @import dplyr
#'
#' @keywords internal

valid_measurement <- function(aggregated_data, status_list) {
  valid_vec <- status_list[["Valid"]]
  valid_thresh <- status_list[["Valid threshold"]]
  ko_vec <- status_list[["KO"]]
  filter_status <- function(vec) {
    status <- sum(vec %in% valid_vec) / length(vec) > valid_thresh
    if (status) {
      if (length(ko_vec) == 0) {
        status <- !any(vec %in% ko_vec)
      }
    }
    return(status)
  }
  aggregated_data <- mutate(aggregated_data,
                            Valid_Replicates = filter_status(Status),
                            .after = Status)
  return(aggregated_data)
}
