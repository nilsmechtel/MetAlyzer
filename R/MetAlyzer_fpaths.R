#' @title Get example data file path
#'
#' @description This function returns the extraction_data.xlsx file path.
#'
#' @return extraction_data.xlsx file path
#' @export
#'
#' @examples
#' fpath <- extraction_data()
extraction_data <- function() {
  system.file("extdata", "extraction_data.xlsx", package = "MetAlyzer")
}


#' @title Get example data file path
#'
#' @description This function returns the treatment_data.xlsx file path.
#'
#' @return treatment_data.xlsx file path
#' @export
#'
#' @examples
#' fpath <- treatment_data()
treatment_data <- function() {
  system.file("extdata", "treatment_data.xlsx", package = "MetAlyzer")
}


#' @title Get polarity file path
#'
#' @description This function returns the polarity.csv file path.
#'
#' @return polarity.csv file path
#' @export
#'
#' @examples
#' fpath <- polarity()
polarity <- function() {
  system.file("extdata", "polarity.csv", package = "MetAlyzer")
}


#' @title Get pathway file path
#'
#' @description This function returns the pathway.xlsx file path.
#'
#' @return pathway.xlsx file path
#' @export
#'
#' @examples
#' fpath <- pathway()
pathway <- function() {
  system.file("extdata", "pathway.xlsx", package = "MetAlyzer")
}

#' @title Get example data xl file path
#'
#' @description This function returns the MetIDQ_XL.xlsx file path.
#'
#' @return MetIDQ_XL.xlsx file path
#' @export
#'
#' @examples
#' fpath <- xl_data()
xl_data <- function() {
  system.file("extdata", "MetIDQ_XL.xlsx", package = "MetAlyzer")
}
