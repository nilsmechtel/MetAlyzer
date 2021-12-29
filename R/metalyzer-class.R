library(xlsx)
library(dplyr)


#' A S4 class to read and analyze MetIDQ output
#'
#' @slot file_path A length-one character vector
#' @slot sheet A length-one numeric vector
#' @slot metabolites A character vector with all 630 measured metabolites and optional 234 additional metabolism indicators
#' @slot raw_data A data frame containing all raw measurements; dimension: # samples x # metabolites
#' @slot quant_status A data frame containing the quantification status of all measurements; dimension: # samples x # metabolites
#' @slot meta_data A data frame containing any meta data; dimension: # samples x # meta variables
#' @slot full_sheet A matrix containing the un-sliced Excel sheet
#' @slot data_ranges A length-six numeric list with rows and columns information for slicing
#'
#' @return
#' @export
#'
#' @examples
#' obj <- new("MetAlyzer")
#' obj <- initialize(obj, file_path="results.xlsx")
#' show(obj)
#' obj <- readData(obj)
#' obj <- readQuantQtatus(obj)
#' obj <- filterClasses(obj)
#' summariseQuantData(obj)

setClass("MetAlyzer",
         slots=list(
           ## input ##
           file_path="character",
           sheet="numeric",
           ## output ##
           metabolites="character",
           raw_data="data.frame",
           quant_status="data.frame",
           meta_data="data.frame",
           ## intern ##
           full_sheet="matrix",
           data_ranges="list"
         )
)


#' Initialize
#'
#' This function initializes file path and sheet number
#' @param object MetAlyzer object
#' @param file_path A length-one character vector with the file path
#' @param sheet A length-one numeric vector with the sheet number
#'
#' @return
#' @export
#'
#' @examples

setGeneric("initialize", function(object, file_path, sheet=1) standardGeneric("initialize"))
setMethod("initialize",
          "MetAlyzer",
          function(object, file_path, sheet=1) {
            object@file_path <- file_path
            object@sheet <- sheet
            return(object)
          })


#' Show a MetAlyzer object
#'
#' This function shows a summary of MetAlyzer slot values
#' @param object MetAlyzer object
#'
#' @return
#' @export
#'
#' @examples

# setGeneric("show", function(object) standardGeneric("show"))
setMethod("show",
          "MetAlyzer",
          function(object) {
            show_obj(object)
          }
)


#' Open file and read data
#'
#' This function opens the given MetIDQ output Excel file and extracts metabolites, raw data and meta data
#' @param object MetAlyzer object
#'
#' @return
#' @export
#'
#' @examples

setGeneric("readData", function(object) standardGeneric("readData"))
setMethod("readData",
          "MetAlyzer",
          function(object) {
            object@full_sheet <- open_file(object)
            object@data_ranges <- get_data_range(object)
            object@metabolites <- read_metabolties(object)
            object@raw_data <- read_raw_data(object)
            object@meta_data <- read_meta_data(object)
            return(object)
          }
)


#' Read quantification status
#'
#' This function gets the background color of each cell of the opened Excel file and assigns it to the corresponding quantification status
#' @param object MetAlyzer object
#'
#' @return
#' @export
#'
#' @examples

setGeneric("readQuantStatus", function(object) standardGeneric("readQuantStatus"))
setMethod("readQuantStatus",
          "MetAlyzer",
          function(object) {
            df_BG <- read_BG_color(object)
            object@quant_status <- assign_quant_status(df_BG)
            return(object)
          }
)


#' Filter metabolites
#'
#' This function filters out certain classes of the metabolites vector and raw_data and quant_status columns
#' @param object MetAlyzer object
#' @param class_name A class to be filtered out
#'
#' @return
#' @export
#'
#' @examples

setGeneric("filterClasses", function(object, class_name="Metabolism Indicators") standardGeneric("filterClasses"))
setMethod("filterClasses",
          "MetAlyzer",
          function(object, class_name) {
            return(filter_classes(object, class_name))
          }
)


#' Summarize quantification status
#'
#' This function lists the proportion of LODs, LOQs, Valids, calibration range passes and NAs
#' @param object MetAlyzer object
#'
#' @return
#' @export
#'
#' @examples

setGeneric("summariseQuantData", function(object) standardGeneric("summariseQuantData"))
setMethod("summariseQuantData",
          "MetAlyzer",
          function(object) {
            sum_quant_data(object)
          }
)
