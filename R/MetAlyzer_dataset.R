#' Open file and read data
#'
#' This function creates a 'MetAlyzer' object, opens the given 'MetIDQ' output Excel
#' sheet and extracts metabolites, raw data, quantification status and meta data.
#' The column "Sample Type" and the row "Class" are used as anchor cells in the
#' Excel sheet and are therefore a requirement.
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
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzerDataset(file_path = fpath)
#' }

MetAlyzer_dataset <- function(file_path, sheet=1) {
  object <- new("MetAlyzer", file_path = file_path, sheet = sheet)
  full_sheet <- open_file(object)
  data_ranges <- get_data_range(full_sheet)

  metabolites <- slice_metabolties(full_sheet, data_ranges)
  object@metabolites <- list("original" = metabolites,
                             "filtered" = metabolites)
  object@meta_data <- slice_meta_data(full_sheet, data_ranges)
  object@raw_data <- slice_raw_data(full_sheet, data_ranges, metabolites)
  object@quant_status <- read_quant_status(file_path, sheet,
                                           nrow(full_sheet), ncol(full_sheet),
                                           data_ranges, metabolites)
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
#' @param full_sheet full_sheet
#'
#' @keywords internal

get_data_range <- function(full_sheet) {
  if (!"Sample Type" %in% full_sheet) {
    stop('The cell "Class" is missing. Metabolites could not be datected...')
  }
  row_class <- which(full_sheet == "Class") %% nrow(full_sheet) # row of header "Class"
  col_class <- ceiling(which(full_sheet == "Class") / nrow(full_sheet)) # column of header "Class"

  if (!"Sample Type" %in% full_sheet) {
    stop('The column "Sample Type" is missing. Data could not be datected...')
  }
  row_sample_type <- which(full_sheet == "Sample Type") %% nrow(full_sheet) # row of "Sample Type"
  col_sample_type <- ceiling(which(full_sheet == "Sample Type") / nrow(full_sheet)) # column of "Sample Type"

  rows_data <- which(full_sheet[, col_sample_type] == "Sample") # rows
  names(rows_data) <- NULL
  cols_data <- (col_class+1):ncol(full_sheet) # columns

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
#' @param full_sheet full_sheet
#' @param data_ranges data_ranges
#'
#' @keywords internal

slice_metabolties <- function(full_sheet, data_ranges) {
  ## metabolites are a row above classes
  metabolites <- full_sheet[data_ranges[["class_row"]]-1,
                            data_ranges[["data_cols"]]]
  classes <- full_sheet[data_ranges[["class_row"]],
                        data_ranges[["data_cols"]]]
  names(metabolites) <- classes
  return(metabolites)
}


#' Slice raw data
#'
#' This function slices measurements from .full_sheet into raw_data.
#'
#' @param full_sheet full_sheet
#' @param data_ranges data_ranges
#' @param metabolites metabolites
#'
#' @keywords internal

slice_raw_data <- function(full_sheet, data_ranges, metabolites) {
  raw_data <- full_sheet[data_ranges[["data_rows"]],
                         data_ranges[["data_cols"]]]
  raw_data <- as.data.frame(raw_data)
  colnames(raw_data) <- metabolites
  raw_data[] <- lapply(raw_data[], as.numeric)
  return(raw_data)
}


#' Slice meta data
#'
#' This function slices meta data from .full_sheet into meta_data.
#'
#' @param full_sheet full_sheet
#' @param data_ranges data_ranges
#'
#' @keywords internal

slice_meta_data <- function(full_sheet, data_ranges) {
  meta_header <- paste(full_sheet[data_ranges[["sample_type_row"]],
                                  1:data_ranges[["class_col"]]])
  meta_data <- as.data.frame(full_sheet[data_ranges[["data_rows"]],
                                        1:data_ranges[["class_col"]]])
  colnames(meta_data) <- meta_header
  meta_data <- meta_data[colSums(is.na(meta_data)) != nrow(meta_data)] # remove full NA columns
  nas <- grep("NA", colnames(meta_data))
  colnames(meta_data)[nas] <- paste0("X", seq_along(nas)) # replace NAs in colnames with X and consecutive numbers
  meta_data$Filter <- TRUE # add extra column for filtering
  return(meta_data)
}



#' Unify hex codes
#'
#' This function gets hex codes of different length and returns them in a
#' unified format.
#'
#' @param hex A 3, 4, 6 or 8 digit hex code
#'
#' @keywords internal

unify_hex <- function(hex) {
  upper_hex <- toupper(hex)
  hex_digits <- stringr::str_extract_all(upper_hex, "[0-9A-F]")[[1]]
  n_digits <- length(hex_digits)
  if (n_digits == 6) {
    hex_code <- paste0("#", paste(hex_digits, collapse = ""))
  } else if (n_digits == 8) {
    hex_code <- paste0("#", paste(hex_digits[3:8], collapse = ""))
  } else if (n_digits == 3) {
    hex_code <- paste0("#", hex_digits[1], hex_digits[1],
                      hex_digits[2], hex_digits[2],
                      hex_digits[3], hex_digits[3])
  } else if (n_digits == 4) {
    hex_code <- paste0("#", hex_digits[2], hex_digits[2],
                      hex_digits[3], hex_digits[3],
                      hex_digits[4], hex_digits[4])
  } else {
    hex_code <- ""
  }
  return(hex_code)
}


#' Read quantification status
#'
#' This function gets the background color of each cell in .full_sheet and
#' assigns it to the corresponding quantification status.
#'
#' @param file_path file_path
#' @param sheet sheet
#' @param n_row n_row
#' @param n_col n_col
#' @param data_ranges data_ranges
#' @param metabolites metabolites
#'
#' @importFrom openxlsx loadWorkbook getSheetNames
#'
#' @keywords internal

read_quant_status <- function(file_path, sheet, n_row, n_col,
                              data_ranges, metabolites) {
  status_list <- list(
    "Valid" = c("#B9DE83", "#00CD66"),
    "LOQ" = c("#B2D1DC","#7FB2C5", "#87CEEB"), # <LLOQ or > ULOQ
    "LOD" = c("#A28BA3", "#6A5ACD"), # < LOD
    "ISTD Out of Range" = c("#FFF099", "#FFFF33"),
    "Invalid" = "#FFFFCC",
    "Incomplete" = c("#CBD2D7", "#FFCCCC")
  )
  wb <- loadWorkbook(file = file_path)
  sheet_name <- getSheetNames(file = file_path)[sheet]
  styles <- wb$styleObjects
  mat_BG <- matrix(NA, ncol = n_col, nrow = n_row)
  for (x in styles) { # iterate through different cell styles in the sheet
    if (x$sheet == sheet_name) {
      rgb <- x$style$fill$fillFg
      rgb <- ifelse(length(rgb) == 0, "", rgb)
      hex_code <- unify_hex(rgb)
      matching_status <- names(status_list)[sapply(status_list, function(colors) {hex_code %in% colors})]
      if (length(matching_status) > 0) {
        if (hex_code != rgb) {
          cat("Converting color code:",
              paste0("\"", rgb, "\""), "->", paste0("\"", hex_code, "\""), "\n")
        }
        status <- matching_status[1] # take the first matching status
        rows <- x$rows # row indices
        cols <- x$cols # column indices
        for (i in 1:length(rows)) {
          mat_BG[rows[i], cols[i]] <- status
        }
      }
    }
  }
  quant_status <- as.data.frame(mat_BG)[data_ranges[["data_rows"]],
                                        data_ranges[["data_cols"]]]
  colnames(quant_status) <- metabolites
  quant_status[] <- lapply(quant_status, function(x) {
    factor(x, levels = names(status_list))
  })
  return(quant_status)
}

