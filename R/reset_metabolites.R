#' Reset metabolites
#'
#' This method resets the filtering of metabolites.
#'
#' @param object MetAlyzer object
#'
#' @keywords internal

reset_metabolites <- function(object) {
  original <- object@metabolites[["original"]]
  filtered <- object@metabolites[["filtered"]]
  diff <- length(original) - length(filtered)
  if (diff > 0) {
    cat(paste("Restoring", diff, "metabolite(s).\n"))
    object@metabolites[["filtered"]] <- original
  }
  return(object)
}
