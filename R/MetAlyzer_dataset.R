#' @title Open file and read data
#'
#' @description This function creates a SummarizedExperiment (SE) from the given
#' 'MetIDQ' output Excel sheet: metabolites (rowData), meta data (colData),
#' concentration data (assay), quantification status(assay)
#' The column "Sample Type" and the row "Class" are used as anchor cells in the
#' Excel sheet and are therefore a requirement.
#'
#' @param file_path A character specifying the file path to the Excel file.
#' @param sheet A numeric index specifying which sheet of the Excel file to use.
#' @param status_list A list of HEX color codes for each quantification status.
#' @param silent If TRUE, mute any print command.
#' @return A Summarized Experiment object
#' @export
#'
#' @examples
#' metalyzer_se <- MetAlyzer_dataset(file_path = example_extraction_data())
MetAlyzer_dataset <- function(
    file_path,
    sheet = 1,
    status_list = list(
      "Valid" = c("#B9DE83", "#00CD66"),
      "LOQ" = c("#B2D1DC", "#7FB2C5", "#87CEEB"), # <LLOQ or > ULOQ
      "LOD" = c("#A28BA3", "#6A5ACD"), # < LOD
      "ISTD Out of Range" = c("#FFF099", "#FFFF33"),
      "Invalid" = "#FFFFCC",
      "Incomplete" = c("#CBD2D7", "#FFCCCC")
    ),
    silent = FALSE) {
  # Print MetAlyzer logo
  if (silent == FALSE) {
    metalyzer_ascii_logo()
  }

  # Open MetIDQ Excel sheet
  starter_list <- list(
    "file_path" = as.character(file_path),
    "sheet" = as.numeric(sheet)
  )
  full_sheet <- open_file(starter_list)
  data_ranges <- get_data_range(full_sheet)

  # Extract metabolites, meta data and concentration values
  metabolites <- slice_metabolites(full_sheet, data_ranges)
  meta_data <- slice_meta_data(full_sheet, data_ranges)
  conc_values <- slice_conc_values(full_sheet, data_ranges, metabolites)

  # Read quantification status
  quant_status <- read_quant_status(
    starter_list = starter_list,
    sheet_dim = c(nrow(full_sheet), ncol(full_sheet)),
    data_ranges = data_ranges,
    metabolites = metabolites,
    status_list = status_list,
    silent = silent
  )

  # Aggregate data and add it to the metadata of SE object
  aggregated_data <- aggregate_data(
    metabolites = metabolites,
    meta_data = meta_data,
    conc_values = conc_values,
    quant_status = quant_status,
    status_vec = names(status_list)
  )

  # Fill Summarized Experiment
  # rowData: features
  # colData: samples
  # assays: conc_values, quant_status
  # metadata: file_path, sheet, aggregated_data

  rowData <- data.frame(
    "metabolic_classes" = names(metabolites),
    row.names = metabolites
  )
  colData <- meta_data
  assays <- list(
    "conc_values" = t(conc_values),
    "quant_status" = t(quant_status)
  )
  metadata <- list(
    "file_path" = file_path,
    "sheet_index" = sheet,
    "status_list" = status_list,
    "aggregated_data" = aggregated_data
  )
  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assays,
    colData = colData,
    rowData = rowData,
    metadata = metadata,
  )

  # Print summary of conc_values and quant_status
  if (silent == FALSE) {
    summarizeConcValues(se)
    summarizeQuantData(se)
  }

  return(se)
}

#' @title MetAlyzer logo
#'
#' @description This function prints "MetAlyzer" as an ASCII logo
#'
#' @keywords internal
metalyzer_ascii_logo <- function() {
  # http://patorjk.com/software/taag/#p=display&f=3D-ASCII&t=MetAlyzer
  # Also try: Rectangles
  line1 <- "\n"
  line2a <- " _____ ______   _______  _________  ________  ___           ___   "
  line2b <- " ___ ________  _______   ________"
  line3a <- "|\\   _ \\  _   \\|\\  ___ \\|\\___   ___\\\\   __  \\|\\  \\     "
  line3b <- "    |\\  \\  /  /|\\_____  \\|\\  ___ \\ |\\   __  \\"
  line4a <- "\\ \\  \\\\\\__\\ \\  \\ \\   __/\\|___ \\  \\_\\ \\  \\|\\  \\ \\"
  line4b <- "  \\        \\ \\  \\/  / /\\|___/  /\\ \\   __/|\\ \\  \\|\\  \\"
  line5a <- " \\ \\  \\\\|__| \\  \\ \\  \\_|/__  \\ \\  \\ \\ \\   __  \\ \\  "
  line5b <- "\\        \\ \\    / /     /  / /\\ \\  \\_|/_\\ \\   _  _\\"
  line6a <- "  \\ \\  \\    \\ \\  \\ \\  \\_|\\ \\  \\ \\  \\ \\ \\  \\ \\  \\"
  line6b <- " \\  \\____    \\/   / /     /  /_/__\\ \\  \\_|\\ \\ \\  \\\\  \\"
  line6c <- "| "
  line7a <- "   \\ \\__\\    \\ \\__\\ \\_______\\  \\ \\__\\ \\ \\__\\ \\__\\ "
  line7b <- "\\_______\\__/   / /     |\\________\\ \\_______\\ \\__\\\\ _\\ "
  line8a <- "    \\|__|     \\|__|\\|_______|   \\|__|  \\|__|\\|__|\\|_______|"
  line8b <- "\\____/ /       \\|_______|\\|_______|\\|__|\\|__|"
  line9 <- "                                                          \\|____|/"
  line10 <- "\n"

  logo_string <- c(
    line1,
    paste(line2a, line2b, sep = ""),
    paste(line3a, line3b, sep = ""),
    paste(line4a, line4b, sep = ""),
    paste(line5a, line5b, sep = ""),
    paste(line6a, line6b, line6c, sep = ""),
    paste(line7a, line7b, sep = ""),
    paste(line8a, line8b, sep = ""),
    line9,
    line10
  )
  cat(logo_string, sep = "\n")
}


#' @title Open Excel file
#'
#' @description This function opens the given MetIDQ output Excel file and reads the full
#' given sheet.
#'
#' @param starter_list contains the file path and the sheet index
#'
#' @keywords internal
open_file <- function(starter_list) {
  full_sheet <- openxlsx::read.xlsx(starter_list$file_path,
    sheet = starter_list$sheet,
    colNames = FALSE,
    skipEmptyRows = FALSE,
    skipEmptyCols = FALSE
  )
  full_sheet <- as.matrix(full_sheet)
  full_sheet[full_sheet == "NA"] <- NA # all entries are strings
  return(full_sheet)
}


#' @title Get data range
#'
#' @description This function extracts rows and column indices to slice the .full_sheet into
#' metabolites, conc_values and meta_data.
#'
#' @param full_sheet full_sheet
#'
#' @keywords internal
get_data_range <- function(full_sheet) {
  if (!"Sample Type" %in% full_sheet) {
    stop('The cell "Class" is missing. Metabolites could not be detected...')
  }
  # row of header "Class"
  row_class <- which(full_sheet == "Class") %% nrow(full_sheet)
  # column of header "Class"
  col_class <- ceiling(which(full_sheet == "Class") / nrow(full_sheet))

  if (!"Sample Type" %in% full_sheet) {
    stop('The column "Sample Type" is missing. Data could not be detected...')
  }
  # row of "Sample Type"
  row_sample_type <- which(full_sheet == "Sample Type") %% nrow(full_sheet)
  # column of "Sample Type"
  col_sample_type <- ceiling(
    which(full_sheet == "Sample Type") / nrow(full_sheet)
  )

  rows_data <- which(full_sheet[, col_sample_type] == "Sample") # rows
  names(rows_data) <- NULL
  cols_data <- (col_class + 1):ncol(full_sheet) # columns

  data_ranges <- list(
    "class_row" = row_class,
    "class_col" = col_class,
    "sample_type_row" = row_sample_type,
    "sample_type_col" = col_sample_type,
    "data_rows" = rows_data,
    "data_cols" = cols_data
  )
  return(data_ranges)
}


#' @title Slice metabolites
#'
#' @description This function extracts metabolites with their corresponding
#' metabolite class from .full_sheet into metabolites.
#'
#' @param full_sheet full_sheet
#' @param data_ranges data_ranges
#'
#' @keywords internal
slice_metabolites <- function(full_sheet, data_ranges) {
  ## metabolites are a row above classes
  metabolites <- full_sheet[
    data_ranges$class_row - 1,
    data_ranges$data_cols
  ]
  classes <- full_sheet[
    data_ranges$class_row,
    data_ranges$data_cols
  ]
  names(metabolites) <- classes
  return(metabolites)
}


#' @title Slice concentration values
#'
#' @description This function slices measurements from .full_sheet into
#' conc_values.
#'
#' @param full_sheet full_sheet
#' @param data_ranges data_ranges
#' @param metabolites metabolites
#'
#' @keywords internal
slice_conc_values <- function(full_sheet, data_ranges, metabolites) {
  conc_values <- full_sheet[
    data_ranges$data_rows,
    data_ranges$data_cols
  ]
  conc_values <- as.data.frame(conc_values)
  colnames(conc_values) <- metabolites
  conc_values[] <- lapply(conc_values[], as.numeric)
  return(conc_values)
}


#' @title Slice meta data
#'
#' @description This function slices meta data from .full_sheet into meta_data.
#'
#' @param full_sheet full_sheet
#' @param data_ranges data_ranges
#'
#' @keywords internal
slice_meta_data <- function(full_sheet, data_ranges) {
  meta_header <- paste(full_sheet[
    data_ranges[["sample_type_row"]],
    1:data_ranges[["class_col"]]
  ])
  meta_data <- as.data.frame(full_sheet[
    data_ranges$data_rows,
    1:data_ranges[["class_col"]]
  ])
  colnames(meta_data) <- meta_header
  meta_data <- meta_data[
    colSums(is.na(meta_data)) != nrow(meta_data)
  ] # remove full NA columns
  nas <- grep("NA", colnames(meta_data))
  colnames(meta_data)[nas] <- paste0(
    "X", seq_along(nas)
  ) # replace NAs in colnames with X and consecutive numbers
  return(meta_data)
}


#' @title Unify hex codes
#'
#' @description This function gets hex codes of different length and returns
#' them in a unified format.
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
    hex_code <- paste0(
      "#", hex_digits[1], hex_digits[1],
      hex_digits[2], hex_digits[2],
      hex_digits[3], hex_digits[3]
    )
  } else if (n_digits == 4) {
    hex_code <- paste0(
      "#", hex_digits[2], hex_digits[2],
      hex_digits[3], hex_digits[3],
      hex_digits[4], hex_digits[4]
    )
  } else {
    hex_code <- ""
  }
  return(hex_code)
}


#' @title Read quantification status
#'
#' @description This function gets the background color of each cell in 
#' full_sheet and assigns it to the corresponding quantification status.
#'
#' @param starter_list starter_list
#' @param sheet_dim sheet_dim
#' @param data_ranges data_ranges
#' @param metabolites metabolites
#' @param status_list status_list
#' @param silent silent
#'
#' @keywords internal
read_quant_status <- function(
    starter_list,
    sheet_dim,
    data_ranges,
    metabolites,
    status_list,
    silent) {
  wb <- openxlsx::loadWorkbook(file = starter_list$file_path)
  sheet_name <- openxlsx::getSheetNames(
    file = starter_list$file_path
  )[starter_list$sheet]
  styles <- wb$styleObjects
  mat_bg <- matrix(NA, nrow = sheet_dim[1], ncol = sheet_dim[2])
  for (x in styles) { # iterate through different cell styles in the sheet
    if (x$sheet == sheet_name) {
      rgb <- x$style$fill$fillFg
      rgb <- ifelse(length(rgb) == 0, "", rgb)
      hex_code <- unify_hex(rgb)
      matching_status <- names(status_list)[
        sapply(status_list, function(colors) {
          hex_code %in% colors
        })
      ]
      if (length(matching_status) > 0) {
        if (hex_code != rgb && silent == FALSE) {
          cat(
            "Info: Reading color code \"", rgb, "\" as \"", hex_code, "\"",
            "\n",
            sep = ""
          )
        }
        status <- matching_status[1] # take the first matching status
        rows <- x$rows # row indices
        cols <- x$cols # column indices
        for (i in seq_along(rows)) {
          mat_bg[rows[i], cols[i]] <- status
        }
      }
    }
  }
  quant_status <- as.data.frame(mat_bg)[
    data_ranges$data_rows,
    data_ranges$data_cols
  ]
  colnames(quant_status) <- metabolites
  return(quant_status)
}

#' @title Aggregate data
#'
#' @description This function reshapes conc_values, quant_status,
#' metatabolites and sample IDs and combines them into a tibble data frame
#' for filtering with dplyr and plotting with 'ggplot2'. "aggregated_data"
#' is grouped by metabolites.
#'
#' @param metabolites metabolites MetAlyzer object
#' @param meta_data Meta_data of the MetAlyzer object
#' @param conc_values conc_values of a MetAlyzer object
#' @param quant_status quant_status of a MetAlyzer object
#' @param status_vec A vector of quantification status
#' @import dplyr
#' @importFrom rlang .data
#'
#' @keywords internal
aggregate_data <- function(
    metabolites,
    meta_data,
    conc_values,
    quant_status,
    status_vec) {
  classes <- names(metabolites)

  comb_data <- dplyr::mutate(
    conc_values,
    ID = rownames(meta_data),
    .before = 1
  )
  gathered_data <- tidyr::gather(comb_data,
    key = "Metabolite",
    value = "Concentration",
    -"ID"
  )
  gathered_status <- tidyr::gather(quant_status,
    key = "Metabolite",
    value = "Status"
  )

  aggregated_data <- gathered_data %>%
    dplyr::group_by(.data$Metabolite) %>%
    dplyr::mutate(
      Class = sapply(.data$Metabolite, function(x) {
        classes[metabolites == x]
      }),
      .after = .data$Metabolite
    )

  aggregated_data$ID <- factor(aggregated_data$ID,
    levels = rownames(meta_data)
  )
  aggregated_data$Metabolite <- factor(aggregated_data$Metabolite,
    levels = unique(metabolites)
  )
  aggregated_data$Class <- factor(aggregated_data$Class,
    levels = unique(classes)
  )
  aggregated_data$Status <- factor(gathered_status$Status,
    levels = status_vec
  )
  aggregated_data <- arrange(aggregated_data, .data$Metabolite)

  return(droplevels(aggregated_data))
}
