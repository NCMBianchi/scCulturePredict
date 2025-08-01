% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evaluation.R
\name{evaluate_predictions}
\alias{evaluate_predictions}
\title{Evaluate prediction results}
\usage{
evaluate_predictions(seurat_object)
}
\arguments{
\item{seurat_object}{A Seurat object containing prediction results in metadata.}
}
\value{
A list containing:
\itemize{
  \item direct_table: Confusion matrix for direct predictions
  \item threshold_table: Confusion matrix for threshold-based predictions
  \item svm_table: Confusion matrix for SVM predictions
  \item direct_accuracy: Accuracy metrics for direct predictions
  \item svm_accuracy: Accuracy metrics for SVM predictions
  \item chi_direct: Chi-squared test results for direct predictions
  \item chi_svm: Chi-squared test results for SVM predictions
}
}
\description{
Evaluates the performance of different prediction methods by computing confusion matrices,
accuracy metrics, and statistical tests.
}
\details{
This function performs a comprehensive evaluation of prediction results by comparing
predicted cell culture media conditions against actual conditions. It supports multiple
prediction methods including direct similarity-based, threshold-based, and SVM predictions.

The evaluation process includes:
\enumerate{
  \item Computing confusion matrices for each prediction method
  \item Calculating accuracy metrics (total correct, percentage) by sample
  \item Performing chi-squared tests to assess statistical significance
  \item Organizing results into a structured list for easy access
}

Confusion matrices show the relationship between actual and predicted conditions,
revealing which conditions are commonly confused with each other. Accuracy metrics
provide an overall assessment of prediction performance, while chi-squared tests
determine whether the predictions are significantly better than random chance.

The function expects prediction results to be stored in the Seurat object metadata
with the following column names:
\itemize{
  \item "predicted_sample_1": Direct similarity-based predictions
  \item "predicted_sample_2": Threshold-based predictions
  \item "classification_pred": SVM-based predictions
}
}
\examples{
evaluation_results <- evaluate_predictions(seurat_object)
}
