% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{create_evaluation_plots}
\alias{create_evaluation_plots}
\title{Create visualization plots}
\usage{
create_evaluation_plots(
  seurat_object,
  results_dir = "./results",
  color_palette = NULL,
  verbose = TRUE
)
}
\arguments{
\item{seurat_object}{A Seurat object containing analysis results.}

\item{results_dir}{Character string specifying the directory to save plots.
Default is "./results".}

\item{color_palette}{Character vector. Colors to use for the plots
(default: NULL, uses default palette).}

\item{verbose}{Logical. Whether to print progress messages (default: TRUE).}
}
\value{
NULL. Plots are saved to the specified directory.
}
\description{
Creates various visualization plots for the analysis results, including UMAP plots
and prediction accuracy plots.
}
\details{
This function automates the creation of standard visualization plots for scCulturePredict
analysis results. It creates and saves a UMAP plot visualizing the experimental
data, which is useful for exploring cell clustering patterns and relationships.

The function performs the following steps:
\enumerate{
  \item Creates a directory to store results if it doesn't already exist
  \item Generates a UMAP plot showing experimental data colored by sample
  \item Saves the plot as a PNG image to the specified directory
  \item Returns NULL after successful execution
}

The UMAP plot visualizes cells in a reduced two-dimensional space, with colors
representing different sample conditions. This visualization helps identify
clusters of cells with similar expression profiles and assess how well different
conditions separate in the reduced dimension space.

This function is particularly useful for quick generation of standard plots
without having to write custom visualization code for common analysis tasks.
}
\examples{
create_evaluation_plots(seurat_object, results_dir = "./results")
}
