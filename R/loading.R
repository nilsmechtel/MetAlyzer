#' Open Excel file
#'
#' This function opens the given MetIDQ output Excel file
#' @param object MetAlyzer object
#'
#' @return
#'
#' @export
#'
#' @examples

open_file <- function(object) {
  full_sheet <- openxlsx::read.xlsx(object@file_path,
                                    sheet = object@sheet,
                                    colNames = FALSE,
                                    skipEmptyRows=FALSE,
                                    skipEmptyCols = FALSE)
  full_sheet <- as.matrix(full_sheet)
  full_sheet[full_sheet == "NA"] <- NA # all entries are strings
  return(full_sheet)
}


#' Get data range
#'
#' This function extracts rows and column indices to slice the .full_sheet into
#' metabolites, raw_data and meta_data
#' @param object MetAlyzer object
#'
#' @return
#'
#' @export
#'
#' @examples

get_data_range <- function(object) {
  row_class <- which(object@.full_sheet == "Class") %% nrow(object@.full_sheet) # row of header "Class"
  col_class <- ceiling(which(object@.full_sheet == "Class") / nrow(object@.full_sheet)) # column of header "Class"

  row_sample_type <- which(object@.full_sheet == "Sample Type") %% nrow(object@.full_sheet) # row of "Sample Type"
  col_sample_type <- ceiling(which(object@.full_sheet == "Sample Type") / nrow(object@.full_sheet)) # column of "Sample Type"

  rows_data <- which(object@.full_sheet[, col_sample_type] == "Sample") # rows
  names(rows_data) <- NULL
  cols_data <- (col_class+1):ncol(object@.full_sheet) # columns

  data_ranges <- list("class_row"=row_class, "class_col"=col_class,
                      "sample_type_row"=row_sample_type, "sample_type_col"=col_sample_type,
                      "data_rows"=rows_data, "data_cols"=cols_data)
  return(data_ranges)
}


#' Get metabolites
#'
#' This function extracts metabolites with their corresponding class from
#' .full_sheet into metabolites
#' @param object MetAlyzer object
#'
#' @return
#'
#' @export
#'
#' @examples

read_metabolties <- function(object) {
  metabolites <- object@.full_sheet[object@.data_ranges[["class_row"]]-1, object@.data_ranges[["data_cols"]]] # metabolites are a row above classes
  classes <- object@.full_sheet[object@.data_ranges[["class_row"]], object@.data_ranges[["data_cols"]]]
  names(metabolites) <- classes
  return(metabolites)
}


#' Get raw data
#'
#' This function slices measurements from .full_sheet into raw_data
#' @param object MetAlyzer object
#'
#' @return
#'
#' @export
#'
#' @examples

read_raw_data <- function(object) {
  raw_data <- as.data.frame(object@.full_sheet[object@.data_ranges[["data_rows"]], object@.data_ranges[["data_cols"]]])
  colnames(raw_data) <- object@metabolites
  raw_data[] <- sapply(raw_data[], as.numeric)
  return(raw_data)
}


#' Get meta data
#'
#' This function slices meta data from .full_sheet into meta_data
#' @param object MetAlyzer object
#'
#' @return
#'
#' @export
#'
#' @examples

read_meta_data <- function(object) {
  meta_header <- paste(object@.full_sheet[object@.data_ranges[["sample_type_row"]], 1:object@.data_ranges[["class_col"]]])
  meta_data <- as.data.frame(object@.full_sheet[object@.data_ranges[["data_rows"]], 1:object@.data_ranges[["class_col"]]])
  colnames(meta_data) <- meta_header
  meta_data <- meta_data[colSums(is.na(meta_data)) != nrow(meta_data)] # remove full NA columns
  nas <- grep("NA", colnames(meta_data))
  colnames(meta_data)[nas] <- paste0("X", seq_along(nas)) # replace NAs in colnames with X and consecutive numbers
  meta_data$Filter <- TRUE # add extra column for filtering
  return(meta_data)
}


#' Get background color
#'
#' This function gets the background color of each cell in .full_sheet
#' @param object MetAlyzer object
#'
#' @return
#'
#' @export
#'
#' @examples

read_BG_color <- function(object) {
  cat("-------------------------------------\n")
  cat("Reading background color...")
  wb <- xlsx::loadWorkbook(object@file_path)
  sheet <- xlsx::getSheets(wb)[[object@sheet]]
  rows <- xlsx::getRows(sheet)
  cells <- xlsx::getCells(rows)
  styles <- sapply(cells, xlsx::getCellStyle) # get style of each cell
  bg <- sapply(styles, function(style) { # get background color of each cell
    fg  <- style$getFillForegroundXSSFColor()
    rgb <- tryCatch(fg$getRgb(), error = function(e) "")
    rgb <- paste(rgb, collapse = "")
    return(rgb)
  })
  row_index <- as.numeric(unlist(lapply(strsplit(names(bg),split = "\\."),
                                        "[", 1))) # row indexes are at first position
  row_index <- row_index - min(row_index) + 1 # start at 1
  col_index <- as.numeric(unlist(lapply(strsplit(names(bg), split = "\\."),
                                        "[", 2))) # column indexes are at second position
  col_index <- col_index - min(col_index) + 1 # start at 1
  mat_BG <- matrix("", ncol = max(col_index), nrow = max(row_index))
  for (i in 1:length(bg)) { # fill background color matrix
    mat_BG[row_index[i], col_index[i]] <- tolower(bg[i])
  }
  df_BG <- as.data.frame(mat_BG)[object@.data_ranges[["data_rows"]], object@.data_ranges[["data_cols"]]]
  colnames(df_BG) <- object@metabolites
  df_BG[is.na(object@raw_data)] <- NA
  cat(" completed!\n")
  return(df_BG)
}


#' Get quantification status
#'
#' This function assigns the background color of each cell in .full_sheet to the
#' corresponding quantification status
#' @param object MetAlyzer object
#'
#' @return
#'
#' @export
#'
#' @examples

assign_quant_status <- function(df_BG) {
  status_list <- list(
    "Valid" = "00cd66",
    "Smaller Zero" = "???",
    "LOD" = "6a5acd", # < LOD
    "LOQ" = "87ceeb", # <LLOQ or > ULOQ
    "No Intercept" = "???",
    "Missing Measurement" = "???",
    "ISTD Out of Range" = "ffff33",
    "Invalid" = "ffffcc",
    "No measurement" = "???"
  )
  for (status in names(status_list)) {
    df_BG[df_BG == status_list[[status]]] <- status
  }
  return(df_BG)
}

