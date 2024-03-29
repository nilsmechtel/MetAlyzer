% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MetAlyzer_dataset.R
\name{MetAlyzer_dataset}
\alias{MetAlyzer_dataset}
\title{Open file and read data}
\usage{
MetAlyzer_dataset(
  file_path,
  sheet = 1,
  status_list = list(Valid = c("#B9DE83", "#00CD66"), LOQ = c("#B2D1DC", "#7FB2C5",
    "#87CEEB"), LOD = c("#A28BA3", "#6A5ACD"), `ISTD Out of Range` = c("#FFF099",
    "#FFFF33"), Invalid = "#FFFFCC", Incomplete = c("#CBD2D7", "#FFCCCC")),
  silent = FALSE
)
}
\arguments{
\item{file_path}{A character specifying the file path to the Excel file.}

\item{sheet}{A numeric index specifying which sheet of the Excel file to use.}

\item{status_list}{A list of HEX color codes for each quantification status.}

\item{silent}{If TRUE, mute any print command.}
}
\value{
A Summarized Experiment object
}
\description{
This function creates a SummarizedExperiment (SE) from the given
'MetIDQ' output Excel sheet: metabolites (rowData), meta data (colData),
concentration data (assay), quantification status(assay)
The column "Sample Type" and the row "Class" are used as anchor cells in the
Excel sheet and are therefore a requirement.
}
\examples{
metalyzer_se <- MetAlyzer_dataset(file_path = example_extraction_data())
}
