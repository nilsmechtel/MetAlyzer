#' Get pathway file path
#'
#' This function returns the pathway.xlsx file path.
#'
#' @return pathway.xlsx file path
#'
#' @export
#'
#' @examples
#' pathway_file()

pathway_file <- function() {
  system.file("extdata", "pathway.xlsx", package = "MetAlyzer")
}
