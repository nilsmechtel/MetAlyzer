#' A S4 class to read and analyze 'MetIDQ' output
#'
#' @slot file_path A length-one character vector giving the file path
#' @slot sheet A length-one numeric vector giving the sheet index
#' @slot metabolites A character vector with all 630 measured metabolites and
#' optional 234 additional metabolism indicators
#' @slot raw_data A data frame containing all raw measurements;
#' dimension: # samples x # metabolites
#' @slot quant_status A data frame containing the quantification status of all
#' measurements; dimension: # samples x # metabolites
#' @slot meta_data A data frame containing any meta data;
#' dimension: # samples x # meta variables
#' @slot plotting_data A tibble data frame containing reshaped information of raw_data,
#' quant_status and meta_data for plotting with ggplot2
#' @slot .full_sheet A matrix containing the un-sliced Excel sheet
#' @slot .data_ranges A length-six numeric list with rows and columns
#' information for slicing
#' @slot .orig_metabolites A character vector storing the unfiltered metabolites
#' vector
#'
#' @export
#'
#' @examples
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#' obj <- filterMetabolites(obj)
#' show(obj)
#' summariseQuantData(obj)
#'
#' obj <- renameMetaData(obj, Method = Group)
#' obj <- filterMetaData(obj, column = Method, keep = 1:6)
#'
#' \donttest{
#' obj <- createPlottingData(obj, Method, Tissue)
#' obj <- imputePlottingData(obj, Method, Metabolite)
#' obj <- transformPlottingData(obj)
#' obj <- performANOVA(obj, categorical = Method)
#' }

MetAlyzer <- setClass("MetAlyzer",
                      slots = list(
                        ## input ##
                        file_path = "character",
                        sheet = "numeric",
                        ## output ##
                        metabolites = "character",
                        raw_data = "data.frame",
                        quant_status = "data.frame",
                        meta_data = "data.frame",
                        plotting_data = "data.frame",
                        ## internal ##
                        .full_sheet = "matrix",
                        .data_ranges = "list",
                        .orig_metabolites = "character"
                      )
)
