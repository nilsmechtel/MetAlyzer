#' Create plotting data
#'
#' This method reshapes raw_data, quant_status and meta_data and combines them
#' together with basic statistics in a tibble data frame for plotting with
#' ggplot2. plotting_data is grouped by metabolites as well as the selection of
#' additional variables. Statistics are then calculated for each group.
#'
#' @param object MetAlyzer object
#' @param ... A selection of columns from meta_data to add to reshaped data frame
#' @param ungrouped A column from meta_data to add to reshaped data frame that
#' will not be used as grouping variables
#' @param ts A numeric vector of thresholds between 0 and 1 for CV categorization
#' @param valid_vec A character vector containing each quantification status that
#' is considered to be a valid measurement
#' @param t A numeric threshold to determine valid measurements
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' \donttest{
#' fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
#' obj <- MetAlyzerDataset(file_path = fpath)
#'
#' obj <- createPlottingData(obj, Tissue, Group,
#' ungrouped = NULL,
#' ts = c(0.1, 0.2, 0.3),
#' valid_vec = c("Valid", "LOQ"), t = 0.5)
#' }


setGeneric("createPlottingData",
           function(object, ...,
                    ungrouped=NULL,
                    ts=c(0.1, 0.2, 0.3),
                    valid_vec=c("Valid", "LOQ"),
                    t=0.5)
             standardGeneric("createPlottingData")
)

#' @describeIn createPlottingData Create plotting data
setMethod("createPlottingData",
          "MetAlyzer",
          function(object, ..., ungrouped, ts, valid_vec, t) {
            create_plotting_data(object, ...,
                                 ungrouped=deparse(substitute(ungrouped)),
                                 ts = ts,
                                 valid_vec = valid_vec,
                                 t = t)
          }
)


#' Impute plotting data
#'
#' This method imputes zero concentration values (Concentration) with the
#' minimal positive value multiplied by i. If all values are zero or NA, they
#' are set to NA. The imputed values are added to plotting_data in an extra
#' column imp_Conc.
#'
#' @param object MetAlyzer object
#' @param ... Variables to group by
#' @param i A numeric value below 1
#' @param imputeNA Logical value whether to impute NA values
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' # To see an example, please check out the vignette.

setGeneric("imputePlottingData",
           function(object, ...,
                    i=0.2,
                    imputeNA=FALSE)
             standardGeneric("imputePlottingData")
)

#' @describeIn imputePlottingData Impute plotting data
setMethod("imputePlottingData",
          "MetAlyzer",
          function(object, ..., i, imputeNA) {
            impute_plotting_data(object, ..., i = i, imputeNA = imputeNA)
          }
)


#' Transform plotting data
#'
#' This method performs a transformation of imputed concentration values
#' (imp_Conc) with a given funciton. NA values are skipped. The transformed
#' values are added to plotting_data in an extra column transf_Conc.
#'
#' @param object MetAlyzer object
#' @param func A function to transform concentration values with
#' number of samples
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' # To see an example, please check out the vignette.

setGeneric("transformPlottingData",
           function(object,
                    func=log2)
             standardGeneric("transformPlottingData")
)

#' @describeIn transformPlottingData Transform plotting data
setMethod("transformPlottingData",
          "MetAlyzer",
          function(object, func) {
            transform_plotting_data(object, func)
          }
)


#' ANOVA
#'
#' This method performs a one-way ANOVA on the grouped plotting_data (the
#' categorical variable is removed from grouping first). For this, the column
#' valid_replicates must have at least one entry that is TRUE in each group.
#' Otherwise, a vector of NA is returned. A Tukey post-hoc test is then used to
#' determine group names, starting with "A" followed by further letters. These
#' group names are added to plotting_data in the column ANOVA_group. Thereby,
#' metabolites can be identified which are significantly higher in one or more
#' of the categorical variable compared to all other for each metabolite.
#'
#' @param object MetAlyzer object
#' @param categorical A column defining the categorical variable
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' # To see an example, please check out the vignette.

setGeneric("performANOVA",
           function(object, categorical)
             standardGeneric("performANOVA")
)

#' @describeIn performANOVA ANOVA
setMethod("performANOVA",
          "MetAlyzer",
          function(object, categorical) {
            perform_ANOVA(object, deparse(substitute(categorical)))
          }
)


#' Update plotting data
#'
#' This method replaces plotting_data with an updated version.
#'
#' @param object MetAlyzer object
#' @param plotting_data An updated plotting_data tibble data frame
#'
#' @return An updated MetAlyzer object
#'
#' @export
#'
#' @examples
#' # To see an example, please check out the vignette.

setGeneric("setPlottingData",
           function(object, plotting_data)
             standardGeneric("setPlottingData")
)

#' @describeIn setPlottingData Update plotting data
setMethod("setPlottingData",
          "MetAlyzer",
          function(object, plotting_data) {
            object@plotting_data <- plotting_data
            return(object)
          }
)
