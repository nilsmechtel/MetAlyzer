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
    sheet = "numeric"
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
  metabolites <- colData(metalyzer)$filtered
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
  meta_data <- as.data.frame(rowData(metalyzer))
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
  conc_values <- assay(metalyzer, "conc_values")
  filtered_sample <- as.data.frame(rowData(metalyzer))$Filter
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
  quant_status <- assay(metalyzer, "quant_status")
  filtered_sample <- as.data.frame(rowData(metalyzer))$Filter
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
  aggregated_data <- assay(metalyzer, "aggregated_data")
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
            if (length(metadata(metalyzer)$file_path) > 0) {
              s_fp <- strsplit(normalizePath(metadata(metalyzer)$file_path), "/")[[1]]
              file <- tail(s_fp, 1)
              path <- paste(s_fp[seq_len(length(s_fp) - 1)], collapse = "/")
            } else {
              file <- "<empty>"
              path <- "<empty>"
            }
            if (length(metadata(metalyzer)$sheet) > 0) {
              sheet <- metadata(metalyzer)$sheet
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
### <--! Does not work rn because of the "MetAlyzer" class -->

#' Summarize concentration values
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
#' na_metabolites <- summarizeConcValues(metalyzer)
#' na_metabolites
summarizeConcValues <- function(metalyzer) {
  conc_values <- getConcValues(metalyzer)
  nas <- sum(is.na(conc_values))
  total <- nrow(conc_values) * ncol(conc_values)
  n_nas <- colSums(is.na(conc_values))
  na_metabolites <- colData(metalyzer)$original[n_nas > 0] #colnames(conc_values)[n_nas > 0] , also possible but maybe conc_values gets cleaned colnames

  cat("-------------------------------------\n")
  cat("Quantiles:\n")
  print(stats::quantile(conc_values, na.rm = TRUE))
  cat(paste0("\nNAs: ", nas, " (", round(nas / total * 100, 2), "%)\n"))
  cat("-------------------------------------\n")

  return(na_metabolites)
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
#' @import SummarizedExperiment
#' @export
#'
#' @examples
#' metalyzer <- MetAlyzer_dataset(file_path = extraction_data())
#' renameMetaData(metalyzer, Method = `Sample Description`)
#' filterMetaData(metalyzer, !is.na(Tissue))
#'
#' aggregated_data <- aggregateData(metalyzer, Tissue, Method)
#' # or
#' result <- aggregateData(metalyzer, Tissue, Method, inplace = FALSE)
#' metalyzer <- result$metalyzer
#' aggregated_data <- result$aggregated_data
aggregateData <- function(metalyzer,
                          ...,
                          status_list = list(
                            "Valid" = c("Valid", "LOQ"),
                            "Valid threshold" = 0.5,
                            "KO" = c()
                          ),
                          inplace = TRUE) {
  env <- parent.frame()  # Get the parent environment
  var_name <- deparse(substitute(metalyzer))  # Get the name of the metalyzer object

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

    assay(metalyzer, "aggregated_data") <- aggregated_data
  } else {
    cat("It seems no data has been loaded.\n")
    cat("Returning empty data.frame to 'aggregated_data' slot.\n")
  }

  if (inplace) {
    # Assign the modified metalyzer object back to the parent environment
    assign(var_name, metalyzer, envir = env)
    return(aggregated_data)
  } else {
    result <- list("metalyzer" = metalyzer, "aggregated_data" = aggregated_data)
    return(result)
  }
}
