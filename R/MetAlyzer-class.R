#' @title MetAlyzer Class
#'
#' @description A S4 class to read and analyze 'MetIDQ' output
#'
#' @slot file_path A length-one character vector giving the file path.
#' @slot sheet A length-one numeric vector giving the sheet index.
#' @slot metabolites A named list of two character vectors, named "original"
#'   and "filtered", containing all 630 measured metabolites and optional
#'   234 additional metabolism indicators.
#' @slot meta_data A data frame containing any meta data.
#'   Dimension: (samples * meta variables).
#' @slot conc_values A data frame containing measured concentration values
#'   for each metabolite of each sample. Dimension: (samples * metabolites).
#' @slot quant_status A data frame containing the quantification status of
#'   each measurement. Dimension: (samples * metabolites).
#' @slot aggregated_data A tibble data frame combining selected variables
#'   from meta_data, metabolite names and classes, conc_values and
#'   quant_status in separate columns.
#'
#' @name MetAlyzer
#' @docType class
#' @export
MetAlyzer <- setClass(
  "MetAlyzer",
  slots = list(
    file_path = "character",
    sheet = "numeric",
    meta_data = "data.frame",
    metabolites = "list",
    conc_values = "data.frame",
    quant_status = "data.frame",
    aggregated_data = "data.frame"
  )
)


# === Getters for MetAlyzer Class Slots ===

#' @title Get metabolites
#'
#' @description This function returns the filtered metabolites vector.
#'
#' @param metalyzer A MetAlyzer object.
#' @return The filtered metabolites vector.
#' @export
#' 
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#' 
#' metabolites <- getMetabolites(metalyzer)
#' metabolites[1:10]
getMetabolites <- function(metalyzer) {
  metabolites <- metalyzer@metabolites[["filtered"]]
  return(metabolites)
}


#' @title Get meta data
#'
#' @description This function returns the meta_data data frame
#' (samples * meta variables) with filtered samples and all columns.
#'
#' @param metalyzer MetAlyzer object
#' @return The filtered meta_data data frame
#' @import dplyr
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' meta_data <- getMetaData(metalyzer)
#' head(meta_data)
getMetaData <- function(metalyzer) {
  meta_data <- metalyzer@meta_data
  if (nrow(meta_data) > 0) {
    meta_data <- meta_data %>%
      dplyr::filter(Filter) %>%
      dplyr::select(-Filter)
  }
  return(meta_data)
}


#' @title Get concentration values
#'
#' @description This function returns the concentration data frame
#' (samples * metabolites) with filtered samples and metabolites.
#'
#' @param metalyzer MetAlyzer object
#' @return The concentration data frame
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' conc_values <- getConcValues(metalyzer)
#' head(conc_values, c(5, 5))
getConcValues <- function(metalyzer) {
  conc_values <- metalyzer@conc_values
  filtered_sample <- metalyzer@meta_data$Filter
  filtered_metabolites <- getMetabolites(metalyzer)
  if (nrow(conc_values) > 0) {
    conc_values <- conc_values[
      filtered_sample,
      filtered_metabolites
    ]
  }
  return(conc_values)
}


#' @title Get quantification status
#'
#' @description This function returns the quant_status data frame (samples *
#'   metabolites) with filtered samples and metabolites.
#'
#' @param metalyzer MetAlyzer object
#' @return The quant_status data frame
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' quant_status <- getQuantStatus(metalyzer)
#' head(quant_status, c(5, 5))
getQuantStatus <- function(metalyzer) {
  quant_status <- metalyzer@quant_status
  filtered_sample <- metalyzer@meta_data$Filter
  filtered_metabolites <- getMetabolites(metalyzer)
  if (nrow(quant_status) > 0) {
    quant_status <- quant_status[
      filtered_sample,
      filtered_metabolites
    ]
  }
  return(quant_status)
}

#' @title Get aggregated data
#'
#' @description This function returns the aggregated data.
#'
#' @param metalyzer A MetAlyzer object.
#' @return The aggregated data.
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#' renameMetaData(metalyzer, Method = `Sample Description`)
#' filterMetaData(metalyzer, !is.na(Tissue))
#' aggregateData(metalyzer, Tissue, Method)
#'
#' aggregated_data <- getAggregatedData(metalyzer)
getAggregatedData <- function(metalyzer) {
  aggregated_data <- metalyzer@aggregated_data
  return(aggregated_data)
}

# === Display MetAlyzer class ===

#' @title Show a 'MetAlyzer' object
#'
#' @description This function shows a summary of the 'MetAlyzer' slot values.
#'
#' @param metalyzer MetAlyzer object
#' @return A summary of the MetAlyzer object
#' @importFrom methods show
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' show(metalyzer)
#' # or
#' metalyzer

setMethod("show",
          "MetAlyzer",
          function(metalyzer) {
            if (length(metalyzer@file_path) > 0) {
              s_fp <- strsplit(normalizePath(metalyzer@file_path), "/")[[1]]
              file <- tail(s_fp, 1)
              path <- paste(s_fp[seq_len(length(s_fp) - 1)], collapse = "/")
            } else {
              file <- "<empty>"
              path <- "<empty>"
            }
            if (length(metalyzer@sheet) > 0) {
              sheet <- metalyzer@sheet
            } else {
              sheet <- "<empty>"
            }
            meta_data <- getMetaData(metalyzer)
            metabolites <- getMetabolites(metalyzer)

            cat("-------------------------------------\n")
            cat("\"MetAlyzer\" object:\n")
            cat("File name: ", file, "\n")
            cat("Sheet: ", sheet, "\n")
            cat("File path: ", path, "\n")
            cat("Metabolites: ", length(metabolites), "\n")
            cat("Classes: ", length(unique(names(metabolites))), "\n")
            if (length(metabolites) > 0) {
              cat("Including metabolism indicators: ",
                  "Metabolism Indicators" %in% names(metabolites), "\n")
            }
            cat("Number of samples: ", nrow(meta_data), "\n")
            if (ncol(meta_data) > 0) {
              cat("Columns meta data: ",
                  paste(colnames(meta_data), collapse = "; "), "\n")
            }
            cat("Aggregated data available: ",
                nrow(getAggregatedData(metalyzer)) > 0,
                "\n"
            )
            cat("-------------------------------------\n")
          })


#' Summarise concentration values
#'
#' This function prints quantiles and NAs of raw data.
#'
#' @param metalyzer MetAlyzer object
#' @return A vector with metabolites that contain NAs
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' na_metabolites <- summariseConcValues(metalyzer)
#' na_metabolites

summariseConcValues <- function(metalyzer) {
  conc_values <- getConcValues(metalyzer)
  nas <- sum(is.na(conc_values))
  total <- nrow(conc_values) * ncol(conc_values)
  n_nas <- colSums(is.na(conc_values))
  na_metabolites <- colnames(conc_values)[n_nas > 0]

  cat("-------------------------------------\n")
  cat("Quantiles:\n")
  print(stats::quantile(conc_values, na.rm = TRUE))
  cat(paste0("\nNAs: ", nas, " (", round(nas / total * 100, 2), "%)\n"))
  cat("-------------------------------------\n")

  return(na_metabolites)
}


#' @title Summarise quantification status
#'
#' @description This function lists the number of each quantification status and
#' its percentage.
#'
#' @param metalyzer MetAlyzer object
#' @return A status_list of metabolite vectors for each quantification status
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' status_list <- summariseQuantData(metalyzer)
#' names(status_list)
summariseQuantData <- function(metalyzer) {
  # Print number of quantification status
  print_number <- function(quant_status, status, total) {
    number <- sum(quant_status == status, na.rm = TRUE)
    cat(paste0(
      status, ": ",
      number, " (", round(number / total * 100, 2), "%)\n"
    ))
  }

  # Get all metabolites that have the quantification status at least once
  status_metabolites <- function(quant_status, status) {
    if (status == "NA") {
      n <- colSums(is.na(quant_status))
    } else {
      n <- colSums(quant_status == status)
    }
    metabolites <- colnames(quant_status)[n > 0]
    return(metabolites)
  }

  quant_status <- getQuantStatus(metalyzer)
  nas <- sum(is.na(quant_status))
  total <- nrow(quant_status) * ncol(quant_status)
  status_vec <- levels(quant_status[, 1])
  status_list <- list()
  cat("-------------------------------------\n")
  for (status in status_vec[which(status_vec %in% unlist(quant_status))]) {
    print_number(quant_status, status, total)
    status_list[[status]] <- status_metabolites(quant_status, status)
  }
  cat(paste0("NAs: ", nas, " (", round(nas / total * 100, 2), "%)\n"))
  status_list[["NA"]] <- status_metabolites(quant_status, "NA")
  cat("-------------------------------------\n")
  return(status_list)
}


# === Handle Meta Data ===

#' @title Filter meta data
#'
#' @description This function updates the "Filter" column in meta_data to
#' filter out samples.
#'
#' @param metalyzer MetAlyzer object
#' @param ... Use ´col_name´ and condition to filter selected variables.
#' @import dplyr
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' filterMetaData(metalyzer, `Sample Description` %in% 1:6)
#' # or
#' filterMetaData(metalyzer, `Sample Description` != 3)
filterMetaData <- function(metalyzer, ...) {
  meta_data <- metalyzer@meta_data
  aggregated_data <- metalyzer@aggregated_data
  orig_len <- sum(meta_data$Filter)

  conditions <- rlang::enquos(...)
  true_rows <- meta_data %>%
    dplyr::mutate(Index = rownames(meta_data)) %>%
    dplyr::filter(!!!rlang::exprs(!!!conditions)) %>%
    dplyr::select(Index) %>%
    unlist()
  meta_data$Filter[!rownames(meta_data) %in% true_rows] <- FALSE
  metalyzer@meta_data <- meta_data

  if (nrow(aggregated_data) > 0) {
      aggregated_data <- aggregated_data %>%
        dplyr::filter(ID %in% rownames(getMetaData(metalyzer))) %>%
        droplevels()
      metalyzer@aggregated_data <- aggregated_data
  }
  diff <- orig_len - sum(meta_data$Filter)
  if (diff == 1) {
    cat("1 sample was excluded!\n")
  } else if (diff > 1) {
    cat(paste(diff, "samples were excluded!\n"))
  } else {
    cat("No samples were excluded!\n")
  }
  metalyzer <<- metalyzer
}


#' @title Reset meta data
#'
#' @description This function resets the filter of meta_data.
#'
#' @param metalyzer MetAlyzer object
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' resetMetaData(metalyzer)
resetMetaData <- function(metalyzer) {
  filter_col <- metalyzer@meta_data$Filter
  if (any(filter_col == FALSE)) {
    cat(paste("Restoring", sum(filter_col == FALSE), "sample(s).\n"))
  }
  metalyzer@meta_data$Filter <- TRUE
  metalyzer <<- metalyzer
}


#' @title Update meta data
#'
#' @description This function adds another column to filtered meta_data.
#'
#' @param metalyzer MetAlyzer object
#' @param ... Use ´new_col_name = new_column´ to rename selected variables
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' updateMetaData(metalyzer, Date = Sys.Date(), Analyzed = TRUE)
updateMetaData <- function(metalyzer, ...) {
  meta_data <- metalyzer@meta_data
  new_cols <- list(...)

  for (col_name in names(new_cols)) {
    new_col <- new_cols[[col_name]]
    if (inherits(new_col, "factor")) {
      levels <- levels(new_col)
    } else {
      levels <- unique(new_col)
    }
    meta_data[, col_name] <- factor(NA, levels = levels)
    meta_data[, col_name][meta_data$Filter == TRUE] <- new_col
  }
  metalyzer@meta_data <- meta_data
  metalyzer <<- metalyzer
}


#' @title Rename meta data
#'
#' @description This function renames a column of meta_data.
#'
#' @param metalyzer MetAlyzer object
#' @param ... Use new_name = old_name to rename selected variables
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' renameMetaData(metalyzer, Method = `Sample Description`)
renameMetaData <- function(metalyzer, ...) {
  metalyzer@meta_data <- dplyr::rename(metalyzer@meta_data, ...)
  metalyzer <<- metalyzer
}


# === Manage Metabolites ===

#' @title Filter metabolites
#'
#' @description This function filters out certain classes or metabolites
#' of the metabolites vector. If aggregated_data is not empty,
#' metabolites and class will also be filtered here.
#'
#' @param metalyzer MetAlyzer object
#' @param drop_metabolites A character vector defining metabolite classes
#' or individual metabolites to be removed
#' @import dplyr
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' drop_metabolites = c("C0", "C2", "C3", "Metabolism Indicators")
#' filterMetabolites(metalyzer, drop_metabolites)
filterMetabolites <- function(metalyzer,
                              drop_metabolites = c("Metabolism Indicators")) {

  metabolites <- metalyzer@metabolites[["filtered"]]
  aggregated_data <- metalyzer@aggregated_data
  orig_len <- length(metabolites)

  rm_metabolites <- drop_metabolites[
    which(drop_metabolites %in% metabolites)
  ]
  rm_classes <- drop_metabolites[
    which(drop_metabolites %in% names(metabolites))
  ]

  # Filter metabolite classes
  if (length(rm_classes) > 0) {
    metabolites <- metabolites[
      -which(names(metabolites) %in% rm_classes)
    ]
    if (nrow(aggregated_data) > 0) {
      aggregated_data <- dplyr::filter(
        aggregated_data,
        !(Class %in% rm_classes)
      ) %>%
      droplevels()
    }
  }
  # Filter individual metabolites
  if (length(rm_metabolites) > 0) {
    metabolites <- metabolites[
      -which(metabolites %in% rm_metabolites)
    ]
    if (nrow(aggregated_data) > 0) {
      aggregated_data <- dplyr::filter(
        aggregated_data,
        !(Metabolite %in% rm_metabolites)
      ) %>%
      droplevels()
    }
  }
  diff <- orig_len - length(metabolites)
  if (diff == 1) {
    cat("1 metabolite was removed!\n")
  } else if (diff > 1) {
    cat(paste(diff, "metabolites were removed!\n"))
  } else {
    cat("No metabolites were removed!\n")
  }
  metalyzer@metabolites[["filtered"]] <- metabolites
  metalyzer@aggregated_data <- aggregated_data
  metalyzer <<- metalyzer
}

#' @title Reset metabolites
#'
#' @description This function resets the filtering of metabolites.
#'
#' @param metalyzer MetAlyzer object
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' resetMetabolites(metalyzer)

resetMetabolites <- function(metalyzer) {
  original <- metalyzer@metabolites[["original"]]
  filtered <- metalyzer@metabolites[["filtered"]]
  diff <- length(original) - length(filtered)
  if (diff > 0) {
    cat(paste("Restoring", diff, "metabolite(s).\n"))
    metalyzer@metabolites[["filtered"]] <- original
  }
  metalyzer <<- metalyzer
}

#' @title Drop invalid metabolites
#'
#' @description This function extracts metabolites with a given percentage
#' of valid replicates from aggregated_data and removes them from
#' metabolites and aggregated data.
#'
#' @param metalyzer MetAlyzer object
#' @param valid_perc A numeric lower threshold between 0 and 1 (t < x) to
#' determine valid metabolites with a given percentage of valid replicates.
#' @import dplyr
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#' renameMetaData(metalyzer, Method = `Sample Description`)
#' filterMetaData(metalyzer, !is.na(Tissue))
#' aggregateData(metalyzer, Tissue, Method)
#'
#' dropInvalidMetabolites(metalyzer)

dropInvalidMetabolites <- function(metalyzer, valid_perc = 0) {
  aggregated_data <- getAggregatedData(metalyzer)
  if (nrow(aggregated_data) == 0) {

  } else {
    cat(paste0("A metabolite counts as invalid, if maximum ",
        round(valid_perc*100, 2),
        "% of measurements are considered as valid for it.\n"))

    df <- aggregated_data %>%
      dplyr::group_by(Metabolite) %>%
      dplyr::summarise(n_valid = sum(Valid_Replicates),
                n_tot = n(),
                perc = n_valid/n_tot) %>%
      dplyr::filter(perc <= valid_perc)
    invalid_metabolites <- as.character(df$Metabolite)
    metalyzer <- filterMetabolites(metalyzer, invalid_metabolites)
  }
  metalyzer <<- metalyzer
}


# === Export data ===


#' @title Export filtered raw data as csv
#'
#' @description This function exports the filtered raw data in the CSV format.
#'
#' @param metalyzer MetAlyzer object
#' @param ... Additional columns from meta_data
#' @param file_path file path
#'
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#'
#' exportConcValues(metalyzer, `Sample Description`, Tissue)


exportConcValues <- function(metalyzer,
                            ...,
                            file_path = "metabolomics_data.csv") {
  meta_data <- getMetaData(metalyzer)
  conc_values <- getConcValues(metalyzer)
  df <- bind_cols(select(meta_data, ...),
                  conc_values)
  cat("Number of samples:", nrow(meta_data), "\n")
  cat("Number of Metabolites:", ncol(conc_values), "\n")
  utils::write.csv(x = df,
                   file = file_path,
                   row.names = FALSE)
}


# === Aggregate data ===

#' @title Aggregate data
#'
#' @description This function reshapes concentration, quant_status
#' and meta_data and combines them in a tibble data frame.
#' "aggregated_data" is grouped by metabolites as well as the selection
#' of additional variables from meta_data. The format is suitable for
#' plotting with ggplot2.
#'
#' @param metalyzer MetAlyzer object
#' @param ... A selection of columns from meta_data to add to aggregated data
#' frame that will be used as extra grouping variables.
#' @param status_list A list with three entries.
#' "Valid" holds the vector with a quantification status that is valid.
#' "Valid threshold" holds the minimal percentage threshold.
#' "KO" holds a character vector with quantification status that directly
#' invalidates a set of repetitions once it is encountered.
#' @return The aggregated data tibble data frame
#' @import dplyr
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#' renameMetaData(metalyzer, Method = `Sample Description`)
#' filterMetaData(metalyzer, !is.na(Tissue))
#'
#' aggregated_data <- aggregateData(metalyzer, Tissue, Method)
aggregateData <- function(metalyzer,
                          ...,
                          status_list = list(
                            "Valid" = c("Valid", "LOQ"),
                            "Valid threshold" = 0.5,
                            "KO" = c()
                          )) {
  meta_data <- getMetaData(metalyzer)
  if (nrow(meta_data) > 0) {
    metabolites <- getMetabolites(metalyzer)
    conc_values <- getConcValues(metalyzer)
    quant_status <- getQuantStatus(metalyzer)
    meta_columns <- select(meta_data, ...)

    cat("Reshape and merge data...  ")
    aggregated_data <- reshape_data(
      metabolites,
      meta_columns,
      conc_values,
      quant_status
    )
    cat("finished!\n")
    grouping_vars <- groups(aggregated_data)

    cat(paste0("Calculate valid replicates (groupwise: ",
              paste(grouping_vars, collapse = " * "), ")...  "))
    aggregated_data <- valid_measurement(
      aggregated_data,
      status_list
    )
    cat("finished!\n")

    metalyzer@aggregated_data <- aggregated_data
  } else {
    cat("It seems no data has been loaded.\n")
    cat("Returning empty data.frame to 'aggregated_data' slot.\n")
  }
  metalyzer <<- metalyzer
  return(aggregated_data)
}
