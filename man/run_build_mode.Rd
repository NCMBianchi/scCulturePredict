% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipeline.R
\name{run_build_mode}
\alias{run_build_mode}
\title{Helper function for build mode}
\usage{
run_build_mode(seurat_object, kegg_file, output_dir, progress, pb, verbose)
}
\value{
A list containing the complete results of BUILD mode analysis:
  \describe{
    \item{seurat_object}{Seurat object with predictions and evaluation metadata}
    \item{models}{List containing trained direct and SVM classification models}
    \item{evaluation_results}{List with accuracy metrics, confusion matrices, and performance statistics}
    \item{plots}{Named list of generated visualization plots}
    \item{fingerprints}{Matrix of calculated transcriptomic fingerprints for each culture condition}
    \item{pathway_activities}{Matrix of pathway activity scores used for training}
  }
}
\description{
Helper function for build mode
}
\examples{
\donttest{
# Example requires prepared data files
# result <- run_build_mode(data_dir = "path/to/data",
#                         experiment_id = "exp1",
#                         output_dir = tempdir())
}
\donttest{
# Example requires reference data and query data
# result <- run_predict_mode(ref_data_dir = "path/to/ref",
#                           query_data_dir = "path/to/query",
#                           output_dir = tempdir())
}
}
\keyword{internal}
