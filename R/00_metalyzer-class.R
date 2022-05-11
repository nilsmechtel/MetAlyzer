#' A S4 class to read and analyze 'MetIDQ' output
#'
#' @slot file_path A length-one character vector giving the file path
#' @slot sheet A length-one numeric vector giving the sheet index
#' @slot metabolites A character vector with all 630 measured metabolites and
#' optional 234 additional metabolism indicators
#' @slot meta_data A data frame containing any meta data;
#' dimension: # samples x # meta variables
#' @slot raw_data A data frame containing all raw measurements;
#' dimension: # samples x # metabolites
#' @slot quant_status A data frame containing the quantification status of all
#' measurements; dimension: # samples x # metabolites
#'
#' @export
#'
#' @examples
#' \dontrun{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' data <- MetAlyzer_dataset(file_path = fpath)
#' }

MetAlyzer <- setClass("MetAlyzer",
                      slots = list(
                        file_path = "character",
                        sheet = "numeric",
                        metabolites = "list",
                        meta_data = "data.frame",
                        raw_data = "data.frame",
                        quant_status = "data.frame"
                      )
)
