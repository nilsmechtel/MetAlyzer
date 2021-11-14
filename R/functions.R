## show ##
showObj <- function(object) {
  s_fp <- strsplit(object@file_path, "/")[[1]]
  file <- tail(s_fp, 1)
  path <- paste(s_fp[1:length(s_fp)-1], collapse = "/")
  cat("-------------------------------------\n")
  cat("File name:", file, "\n")
  cat("Sheet:", object@sheet, "\n")
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

# read the full Excel sheet ##
openFile <- function(object) {
  full_sheet <- openxlsx::read.xlsx(object@file_path,
                                    sheet = object@sheet,
                                    colNames = FALSE,
                                    skipEmptyRows=FALSE,
                                    skipEmptyCols = FALSE)
  full_sheet <- as.matrix(full_sheet)
  full_sheet[full_sheet == "NA"] <- NA # all entries are strings
  return(full_sheet)
}

## data range ##
getDataRange <- function(object) {
  row_class <- which(object@full_sheet == "Class") %% nrow(object@full_sheet) # row of header "Class"
  col_class <- ceiling(which(object@full_sheet == "Class") / nrow(object@full_sheet)) # column of header "Class"

  row_sample_type <- which(object@full_sheet == "Sample Type") %% nrow(object@full_sheet) # row of "Sample Type"
  col_sample_type <- ceiling(which(object@full_sheet == "Sample Type") / nrow(object@full_sheet)) # column of "Sample Type"

  rows_data <- which(object@full_sheet[, col_sample_type] == "Sample") # rows
  names(rows_data) <- NULL
  cols_data <- (col_class+1):ncol(object@full_sheet) # columns

  data_ranges <- list("class_row"=row_class, "class_col"=col_class,
                      "sample_type_row"=row_sample_type, "sample_type_col"=col_sample_type,
                      "data_rows"=rows_data, "data_cols"=cols_data)
  return(data_ranges)
}

## metabolites with classes ##
readMetabolties <- function(object) {
  metabolites <- object@full_sheet[object@data_ranges[["class_row"]]-1, object@data_ranges[["data_cols"]]] # metabolites are a row above classes
  classes <- object@full_sheet[object@data_ranges[["class_row"]], object@data_ranges[["data_cols"]]]
  names(metabolites) <- classes
  return(metabolites)
}

## raw data ##
readRawData <- function(object) {
  raw_data <- as.data.frame(object@full_sheet[object@data_ranges[["data_rows"]], object@data_ranges[["data_cols"]]])
  colnames(raw_data) <- object@metabolites
  raw_data[] <- sapply(raw_data[], as.numeric)
  return(raw_data)
}

## meta data ##
readMetaData <- function(object) {
  meta_header <- paste(object@full_sheet[object@data_ranges[["sample_type_row"]], 1:object@data_ranges[["class_col"]]])
  meta_data <- as.data.frame(object@full_sheet[object@data_ranges[["data_rows"]], 1:object@data_ranges[["class_col"]]])
  colnames(meta_data) <- meta_header
  meta_data <- meta_data[colSums(is.na(meta_data)) != nrow(meta_data)] # remove full NA columns
  nas <- grep("NA", colnames(meta_data))
  colnames(meta_data)[nas] <- paste0("X", seq_along(nas)) # replace NAs in colnames with X and consecutive numbers
  return(meta_data)
}

## read background color ##
readBGColor <- function(object) {
  wb <- loadWorkbook(object@file_path)
  sheet1 <- getSheets(wb)[[object@sheet]]
  rows <- getRows(sheet1)
  cells <- getCells(rows)
  styles <- sapply(cells, getCellStyle) # get style of each cell
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
    mat_BG[row_index[i], col_index[i]] <- bg[i]
  }
  df_BG <- as.data.frame(mat_BG)[object@data_ranges[["data_rows"]], object@data_ranges[["data_cols"]]]
  colnames(df_BG) <- object@metabolites
  df_BG[is.na(object@raw_data)] <- NA
  return(df_BG)
}

## assign back ground color to quantification status ##
assignQuantStatus <- function(df_BG) {
  df_BG[df_BG == "6a5acd"] <- "LOD"
  df_BG[df_BG == "87ceeb"] <- "LOQ"
  df_BG[df_BG == "00cd66"] <- "Valid"
  df_BG[df_BG == "ffff33"] <- "Out of calibration range"
  return(df_BG)
}

## filter matabolism indicators ##
filterClasses_f <- function(object, class_name="Metabolism Indicators") {
  if (class_name %in% names(object@metabolites)) {
    not_indicators <- which(names(object@metabolites) != class_name)
    metabolites <- object@metabolites[not_indicators]

    object@raw_data <- object@raw_data[, metabolites]
    object@quant_status <- object@quant_status[, metabolites]
    object@metabolites <- metabolites
  } else {
    print("-------------------------------------")
    print(paste("No", class_name, "to filter! Returning original object."))
  }
  return(object)
}

## summarise qunatification data ##
sumQunatData <- function(object) {
  sum_LOD <- sum(obj@quant_status == "LOD", na.rm = TRUE)
  sum_LOQ <- sum(obj@quant_status == "LOQ", na.rm = TRUE)
  sum_valid <- sum(obj@quant_status == "Valid", na.rm = TRUE)
  sum_val_range <- sum(obj@quant_status == "Out of calibration range", na.rm = TRUE)
  nas <- sum(is.na(obj@quant_status))
  total <- nrow(obj@quant_status)*ncol(obj@quant_status)
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
