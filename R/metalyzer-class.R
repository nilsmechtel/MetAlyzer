library(xlsx)
library(dplyr)

setClass("metalyzer",
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

## show ##
setGeneric("show", function(object) standardGeneric("show"))
setMethod("show",
          "metalyzer",
          function(object) {
            showObj(object)
          }
)

## open file and read data ##
setGeneric("readData", function(object) standardGeneric("readData"))
setMethod("readData",
          "metalyzer",
          function(object) {
            object@full_sheet <- openFile(object)
            object@data_ranges <- getDataRange(object)
            object@metabolites <- readMetabolties(object)
            object@raw_data <- readRawData(object)
            object@meta_data <- readMetaData(object)
            return(object)
          }
)

## read quantification status ##
setGeneric("readQuantQtatus", function(object) standardGeneric("readQuantQtatus"))
setMethod("readQuantQtatus",
          "metalyzer",
          function(object) {
            df_BG <- readBGColor(object)
            object@quant_status <- assignQuantStatus(df_BG)
            return(object)
          }
)

## filter metabolism indicators ##
setGeneric("filterClasses", function(object, class_name="Metabolism Indicators") standardGeneric("filterClasses"))
setMethod("filterClasses",
          "metalyzer",
          function(object, class_name) {
            return(filterClasses_f(object, class_name))
          }
)


## summarise qunatification data ##
setGeneric("summariseQunatData", function(object) standardGeneric("summariseQunatData"))
setMethod("summariseQunatData",
          "metalyzer",
          function(object) {
            sumQunatData(object)
          }
)
