% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/metalyzer.R
\name{summarizeQuantData}
\alias{summarizeQuantData}
\alias{summarizeQuantData,MetAlyzer-method}
\title{Summarize quantification status}
\usage{
summarizeQuantData(object)

\S4method{summarizeQuantData}{MetAlyzer}(object)
}
\arguments{
\item{object}{MetAlyzer object}
}
\value{
A summary of the quantification status
}
\description{
This method lists the number of each quantification status and its
percentage.
}
\section{Methods (by class)}{
\itemize{
\item \code{summarizeQuantData(MetAlyzer)}: Summarize quantification status

}}
\examples{
\dontrun{
fpath <- system.file("extdata", "example_data.xlsx", package = "MetAlyzer")
data <- MetAlyzerDataset(file_path = fpath)

summarizeQuantData(data)
}
}