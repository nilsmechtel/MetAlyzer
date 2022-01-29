#' Open file and read data
#'
#' This function creates a 'MetAlyzer' object, opens the given 'MetIDQ' output Excel
#' sheet and extracts metabolites, raw data, quantification status and meta data.
#'
#' @param file_path file path
#' @param sheet sheet index
#'
#' @return An MetAlyzer object
#'
#' @importFrom methods new
#' @export
#'
#' @examples
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)

MetAlyzerDataset <- function(file_path, sheet=1) {
  object <- new("MetAlyzer", file_path = file_path, sheet = sheet)
  object@.full_sheet <- open_file(object)
  object@.data_ranges <- get_data_range(object)
  object@.orig_metabolites <- slice_metabolties(object)
  object@metabolites <- object@.orig_metabolites
  object@meta_data <- slice_meta_data(object)
  object@raw_data <- slice_raw_data(object)
  object@quant_status <- read_quant_status(object)
  return(object)
}

#' Open Excel file
#'
#' This function opens the given MetIDQ output Excel file and reads the full
#' given sheet.
#'
#' @param object MetAlyzer object
#'
#' @importFrom openxlsx read.xlsx
#'
#' @keywords internal

open_file <- function(object) {
  full_sheet <- read.xlsx(object@file_path,
                          sheet = object@sheet,
                          colNames = FALSE,
                          skipEmptyRows = FALSE,
                          skipEmptyCols = FALSE)
  full_sheet <- as.matrix(full_sheet)
  full_sheet[full_sheet == "NA"] <- NA # all entries are strings
  return(full_sheet)
}


#' Get data range
#'
#' This function extracts rows and column indices to slice the .full_sheet into
#' metabolites, raw_data and meta_data.
#'
#' @param object MetAlyzer object
#'
#' @keywords internal

get_data_range <- function(object) {
  row_class <- which(object@.full_sheet == "Class") %% nrow(object@.full_sheet) # row of header "Class"
  col_class <- ceiling(which(object@.full_sheet == "Class") / nrow(object@.full_sheet)) # column of header "Class"

  row_sample_type <- which(object@.full_sheet == "Sample Type") %% nrow(object@.full_sheet) # row of "Sample Type"
  col_sample_type <- ceiling(which(object@.full_sheet == "Sample Type") / nrow(object@.full_sheet)) # column of "Sample Type"

  rows_data <- which(object@.full_sheet[, col_sample_type] == "Sample") # rows
  names(rows_data) <- NULL
  cols_data <- (col_class+1):ncol(object@.full_sheet) # columns

  data_ranges <- list("class_row" = row_class,
                      "class_col" = col_class,
                      "sample_type_row" = row_sample_type,
                      "sample_type_col" = col_sample_type,
                      "data_rows" = rows_data,
                      "data_cols" = cols_data)
  return(data_ranges)
}


#' Slice metabolites
#'
#' This function extracts metabolites with their corresponding metabolite class
#' from .full_sheet into metabolites.
#'
#' @param object MetAlyzer object
#'
#' @keywords internal

slice_metabolties <- function(object) {
  metabolites <- object@.full_sheet[object@.data_ranges[["class_row"]]-1, object@.data_ranges[["data_cols"]]] # metabolites are a row above classes
  classes <- object@.full_sheet[object@.data_ranges[["class_row"]], object@.data_ranges[["data_cols"]]]
  names(metabolites) <- classes
  return(metabolites)
}


#' Slice raw data
#'
#' This function slices measurements from .full_sheet into raw_data.
#'
#' @param object MetAlyzer object
#'
#' @keywords internal

slice_raw_data <- function(object) {
  raw_data <- object@.full_sheet[object@.data_ranges[["data_rows"]],
                                 object@.data_ranges[["data_cols"]]]
  raw_data <- as.data.frame(raw_data)
  colnames(raw_data) <- object@metabolites
  raw_data[] <- lapply(raw_data[], as.numeric)
  return(raw_data)
}


#' Slice meta data
#'
#' This function slices meta data from .full_sheet into meta_data.
#'
#' @param object MetAlyzer object
#'
#' @keywords internal

slice_meta_data <- function(object) {
  meta_header <- paste(object@.full_sheet[object@.data_ranges[["sample_type_row"]], 1:object@.data_ranges[["class_col"]]])
  meta_data <- as.data.frame(object@.full_sheet[object@.data_ranges[["data_rows"]], 1:object@.data_ranges[["class_col"]]])
  colnames(meta_data) <- meta_header
  meta_data <- meta_data[colSums(is.na(meta_data)) != nrow(meta_data)] # remove full NA columns
  nas <- grep("NA", colnames(meta_data))
  colnames(meta_data)[nas] <- paste0("X", seq_along(nas)) # replace NAs in colnames with X and consecutive numbers
  meta_data$Filter <- TRUE # add extra column for filtering
  return(meta_data)
}


#' Read quantification status
#'
#' This function gets the background color of each cell in .full_sheet and
#' assigns it to the corresponding quantification status.
#'
#' @param object MetAlyzer object
#'
#' @importFrom openxlsx loadWorkbook getSheetNames
#'
#' @keywords internal

read_quant_status <- function(object) {
  status_list <- list(
    "Valid" = "00CD66",
    "LOQ" = "87CEEB", # <LLOQ or > ULOQ
    "LOD" = "6A5ACD", # < LOD
    "ISTD Out of Range" = "FFFF33",
    "Invalid" = "FFFFCC",
    "Incomplete" = "FFCCCC"
  )
  wb <- loadWorkbook(file = object@file_path)
  sheet_name <- getSheetNames(file = object@file_path)[object@sheet]
  styles <- wb$styleObjects
  mat_BG <- matrix(NA, ncol = ncol(object@.full_sheet),
                   nrow = nrow(object@.full_sheet))
  for (x in styles) { # iterate through different styles in the sheet
    if (x$sheet == sheet_name) {
      rgb <- toupper(substr(x$style$fill$fillFg, 3,8))
      rgb <- ifelse(length(rgb) == 0, "", rgb)
      if (rgb %in% status_list) {
        status <- names(which(status_list == rgb))
        rows <- x$rows # row indices
        cols <- x$cols # colum indices
        for (i in 1:length(rows)) {
          mat_BG[rows[i], cols[i]] <- status
        }
      }
    }
  }
  quant_status <- as.data.frame(mat_BG)[object@.data_ranges[["data_rows"]],
                                 object@.data_ranges[["data_cols"]]]
  colnames(quant_status) <- object@metabolites
  # quant_status[is.na(object@raw_data)] <- NA
  quant_status[] <- lapply(quant_status, function(x) {
    factor(x, levels = names(status_list))
  })
  return(quant_status)
}
