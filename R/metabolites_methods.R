#' Filter metabolites
#'
#' This method filters out certain classes or metabolites of the metabolites
#' vector.
#'
#' Note: If both "metabo_vec" and "class_name" arguments are used "metabo_vec"
#' overwrites the "class_name" argument!
#'
#' @param object MetAlyzer object
#' @param class_name A character value defining the class to be removed
#' @param metabo_vec A character vector defining metabolites to be removed
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' obj <- filterMetabolites(obj, class_name = "Metabolism Indicators")
#' # or
#' obj <- filterMetabolites(obj, metabo_vec = c("C0", "C2", "C3"))

setGeneric("filterMetabolites",
           function(object,
                    class_name="Metabolism Indicators",
                    metabo_vec=NULL)
             standardGeneric("filterMetabolites")
)

#' @describeIn filterMetabolites Filter metabolites
setMethod("filterMetabolites",
          "MetAlyzer",
          function(object, class_name, metabo_vec) {
            filter_metabolites(object, class_name, metabo_vec)
          }
)


#' Reset metabolites
#'
#' This method resets the filtering of metabolites.
#'
#' @param object MetAlyzer object
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' obj <- resetMetabolites(obj)

setGeneric("resetMetabolites",
           function(object)
             standardGeneric("resetMetabolites")
)

#' @describeIn resetMetabolites Reset metabolites
setMethod("resetMetabolites",
          "MetAlyzer",
          function(object) {
            reset_metabolites(object)
          }
)

