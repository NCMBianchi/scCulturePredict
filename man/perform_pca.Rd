% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dimensionality_reduction.R
\name{perform_pca}
\alias{perform_pca}
\title{Perform PCA analysis}
\usage{
perform_pca(seurat_object, n_pcs = 50, features = NULL)
}
\arguments{
\item{seurat_object}{A Seurat object containing preprocessed single-cell data.}

\item{n_pcs}{Integer specifying the number of principal components to compute.
Default is 50.}

\item{features}{Character vector specifying which features to use for PCA.
If NULL, uses variable features. Default is NULL.}
}
\value{
A list containing:
\itemize{
  \item seurat_object: The Seurat object with PCA results
  \item variance_explained: Data frame with variance explained by each PC
  \item elbow_plot: ggplot object showing the elbow plot
}
}
\description{
Performs Principal Component Analysis (PCA) on a Seurat object and returns
the results along with variance explained by each principal component.
}
\examples{
pca_results <- perform_pca(seurat_object, n_pcs = 50)
}
