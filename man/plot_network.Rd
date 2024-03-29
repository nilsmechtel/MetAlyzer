% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_network.R
\name{plot_network}
\alias{plot_network}
\title{Plot Pathway Network}
\usage{
plot_network(
  metalyzer_se,
  q_value = 0.05,
  metabolite_text_size = 3,
  connection_width = 0.75,
  pathway_text_size = 6,
  pathway_width = 4,
  scale_colors = c("green", "black", "magenta")
)
}
\arguments{
\item{metalyzer_se}{A Metalyzer object}

\item{q_value}{The q-value threshold for significance}

\item{metabolite_text_size}{The text size of metabolite labels}

\item{connection_width}{The line width of connections between metabolites}

\item{pathway_text_size}{The text size of pathway annotations}

\item{pathway_width}{The line width of pathway-specific connection coloring}

\item{scale_colors}{A vector of length 3 with colors for low, mid and high 
of the gradient.}
}
\value{
ggplot object
}
\description{
This function plots the log2 fold change for each metabolite and visualizes it, in a pathway network.
}
\examples{
metalyzer_se <- MetAlyzer_dataset(file_path = example_mutation_data_xl())
metalyzer_se <- filterMetabolites(
  metalyzer_se,
  drop_metabolites = "Metabolism Indicators"
)
metalyzer_se <- renameMetaData(
  metalyzer_se,
  Mutant_Control = "Sample Description"
)

metalyzer_se <- calculate_log2FC(
  metalyzer_se,
  categorical = "Mutant_Control",
  impute_perc_of_min = 0.2,
  impute_NA = FALSE
)

network <- plot_network(metalyzer_se, q_value = 0.05)
}
