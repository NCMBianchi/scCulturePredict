% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotting.R
\name{plot_scCulture}
\alias{plot_scCulture}
\title{Plot scCulture Results}
\usage{
plot_scCulture(
  scCulture_results,
  plot_type = NULL,
  point_size = 1.5,
  point_alpha = NULL,
  return_data = FALSE
)
}
\arguments{
\item{scCulture_results}{List. Results object returned by scCulture() function.
Must contain seurat_object and mode information.}

\item{plot_type}{Character. Type of plot to generate. For BUILD mode: "accuracy" (default).
For PREDICT mode: "both" (default), "predictions", or "confidence".}

\item{point_size}{Numeric. Size of points in UMAP plots (default: 1.5).}

\item{point_alpha}{Numeric. Transparency of points (default: 0.7 for BUILD, 0.8 for PREDICT).}

\item{return_data}{Logical. If TRUE, returns the data used for plotting instead of plots (default: FALSE).}
}
\value{
For BUILD mode: A ggplot object showing prediction accuracy on UMAP.
For PREDICT mode: A list of ggplot objects containing both predictions and confidence plots (default),
  or a single ggplot object when plot_type is "predictions" or "confidence".
If return_data = TRUE, returns the data frame used for plotting.
}
\description{
Automatically generates appropriate visualizations for scCulture analysis results.
Detects whether results are from BUILD or PREDICT mode and creates corresponding plots.
}
\details{
This function automatically detects the analysis mode from the scCulture results and generates
appropriate visualizations:

**BUILD Mode Plots:**
\itemize{
  \item **Accuracy Plot**: UMAP colored by prediction accuracy (blue = correct, red = incorrect)
  \item Shows how well the trained model performs on the training data
  \item Helps identify regions where the model struggles
}

**PREDICT Mode Plots:**
\itemize{
  \item **Both Plots (default)**: Returns both predictions and confidence plots as a list
  \item **Predictions Plot**: UMAP colored by predicted culture medium (single plot)
  \item **Confidence Plot**: UMAP colored by prediction confidence scores (single plot)
}

The function extracts UMAP coordinates and relevant metadata from the Seurat object
within the scCulture results. It handles missing data gracefully and provides
informative error messages for common issues.
}
\examples{
# Example with mock scCulture results
# Create mock result structure
mock_results <- list(
    predictions = data.frame(
        cell_id = paste0("Cell", seq_len(100)),
        predicted_culture = sample(c("Culture1", "Culture2", "Culture3"), 100, replace = TRUE),
        confidence = runif(100, 0.6, 0.95)
    ),
    umap_coords = data.frame(
        UMAP_1 = rnorm(100),
        UMAP_2 = rnorm(100)
    )
)
class(mock_results) <- "scCulture"

# Create visualization
plot <- plot_scCulture(mock_results, plot_type = "confidence")
}
\seealso{
\code{\link{scCulture}} for generating the results object
}
