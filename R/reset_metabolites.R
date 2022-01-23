#' Reset metabolites
#'
#' This method resets the filtering of metabolites.
#'
#' @param object MetAlyzer object
#'
#' @keywords internal

reset_metabolites <- function(object) {
  orig <- object@.orig_metabolites
  diff <- length(orig) - length(object@metabolites)
  if (diff > 0) {
    cat(paste("Restoring", diff, "metabolite(s).\n"))
  }
  object@metabolites <- orig
  return(object)
}
