#' Evaluate prediction results
#'
#' @description
#' Evaluates the performance of different prediction methods by computing confusion matrices,
#' accuracy metrics, and statistical tests.
#'
#' @details
#' This function performs a comprehensive evaluation of prediction results by comparing
#' predicted cell culture media conditions against actual conditions. It supports multiple
#' prediction methods including direct similarity-based, threshold-based, and SVM predictions.
#'
#' The evaluation process includes:
#' \enumerate{
#'   \item Computing confusion matrices for each prediction method
#'   \item Calculating accuracy metrics (total correct, percentage) by sample
#'   \item Performing chi-squared tests to assess statistical significance
#'   \item Organizing results into a structured list for easy access
#' }
#'
#' Confusion matrices show the relationship between actual and predicted conditions,
#' revealing which conditions are commonly confused with each other. Accuracy metrics
#' provide an overall assessment of prediction performance, while chi-squared tests
#' determine whether the predictions are significantly better than random chance.
#'
#' The function expects prediction results to be stored in the Seurat object metadata
#' with the following column names:
#' \itemize{
#'   \item "predicted_sample_1": Direct similarity-based predictions
#'   \item "predicted_sample_2": Threshold-based predictions
#'   \item "classification_pred": SVM-based predictions
#' }
#'
#' @param seurat_object A Seurat object containing prediction results in metadata.
#'
#' @return A list containing:
#' \itemize{
#'   \item direct_table: Confusion matrix for direct predictions
#'   \item threshold_table: Confusion matrix for threshold-based predictions
#'   \item svm_table: Confusion matrix for SVM predictions
#'   \item direct_accuracy: Accuracy metrics for direct predictions
#'   \item svm_accuracy: Accuracy metrics for SVM predictions
#'   \item chi_direct: Chi-squared test results for direct predictions
#'   \item chi_svm: Chi-squared test results for SVM predictions
#' }
#' @export
#'
#' @examples
#' # Create mock Seurat object with predictions
#' library(Seurat)
#' counts <- matrix(rpois(500, 5), nrow = 50)
#' rownames(counts) <- paste0("Gene", seq_len(50))
#' colnames(counts) <- paste0("Cell", seq_len(10))
#'
#' # Create metadata with actual and predicted samples
#' metadata <- data.frame(
#'   row.names = colnames(counts),
#'   sample = rep(c("A", "B"), each = 5),
#'   predicted_sample_1 = c("A", "A", "B", "A", "A", "B", "B", "B", "A", "B"),
#'   predicted_sample_2 = c("A", "A", "A", "A", "A", "B", "B", "B", "B", "B"),
#'   classification_pred = c("A", "A", "B", "A", "A", "B", "B", "B", "B", "B")
#' )
#' seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
#'
#' # Evaluate predictions
#' results <- evaluate_predictions(seurat_obj)
#'
#' # Check results
#' print(results$direct_accuracy)
evaluate_predictions <- function(seurat_object) {
  # Confusion matrices
  direct_table <- table(seurat_object$sample, seurat_object$predicted_sample_1)
  threshold_table <- table(seurat_object$sample, seurat_object$predicted_sample_2)
  svm_table <- table(seurat_object$sample, seurat_object$classification_pred)

  # Calculate accuracy percentages
  meta_data <- seurat_object@meta.data

  # Direct accuracy calculation
  direct_accuracy_list <- by(meta_data, meta_data$sample, function(x) {
    correct_rate <- mean(x$sample == x$predicted_sample_1, na.rm = TRUE)
    data.frame(
      sample = x$sample[1],
      correct = correct_rate,
      count = nrow(x),
      percent = sprintf("%.1f%%", 100 * correct_rate),
      stringsAsFactors = FALSE
    )
  })
  direct_accuracy <- do.call(rbind, direct_accuracy_list)
  rownames(direct_accuracy) <- NULL

  # SVM accuracy calculation
  svm_accuracy_list <- by(meta_data, meta_data$sample, function(x) {
    correct_rate <- mean(x$sample == x$classification_pred, na.rm = TRUE)
    data.frame(
      sample = x$sample[1],
      correct = correct_rate,
      count = nrow(x),
      percent = sprintf("%.1f%%", 100 * correct_rate),
      stringsAsFactors = FALSE
    )
  })
  svm_accuracy <- do.call(rbind, svm_accuracy_list)
  rownames(svm_accuracy) <- NULL

  # Chi-squared tests
  chi_direct <- tryCatch(chisq.test(direct_table), error = function(e) NULL)
  chi_svm <- tryCatch(chisq.test(svm_table), error = function(e) NULL)

  return(list(
    direct_table = direct_table,
    threshold_table = threshold_table,
    svm_table = svm_table,
    direct_accuracy = direct_accuracy,
    svm_accuracy = svm_accuracy,
    chi_direct = chi_direct,
    chi_svm = chi_svm
  ))
}

#' Create visualization plots
#'
#' @description
#' Creates various visualization plots for the analysis results, including UMAP plots
#' and prediction accuracy plots.
#'
#' @details
#' This function automates the creation of standard visualization plots for scCulturePredict
#' analysis results. It creates and saves a UMAP plot visualizing the experimental
#' data, which is useful for exploring cell clustering patterns and relationships.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates a directory to store results if it doesn't already exist
#'   \item Generates a UMAP plot showing experimental data colored by sample
#'   \item Saves the plot as a PNG image to the specified directory
#'   \item Returns NULL after successful execution
#' }
#'
#' The UMAP plot visualizes cells in a reduced two-dimensional space, with colors
#' representing different sample conditions. This visualization helps identify
#' clusters of cells with similar expression profiles and assess how well different
#' conditions separate in the reduced dimension space.
#'
#' This function is particularly useful for quick generation of standard plots
#' without having to write custom visualization code for common analysis tasks.
#'
#' @param seurat_object A Seurat object containing analysis results.
#' @param results_dir Character string specifying the directory to save plots.
#' Default is "./results".
#' @param color_palette Character vector. Colors to use for the plots
#'   (default: NULL, uses default palette).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return NULL. Plots are saved to the specified directory.
#' @export
#'
#' @examples
#' # Create mock Seurat object with UMAP coordinates
#' library(Seurat)
#' counts <- matrix(rpois(500, 5), nrow = 50)
#' rownames(counts) <- paste0("Gene", seq_len(50))
#' colnames(counts) <- paste0("Cell", seq_len(10))
#'
#' # Create metadata with UMAP coordinates
#' metadata <- data.frame(
#'   row.names = colnames(counts),
#'   UMAP_1 = rnorm(10),
#'   UMAP_2 = rnorm(10),
#'   sample = rep(c("A", "B"), each = 5)
#' )
#' seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
#'
#' # Create plots in temporary directory
#' temp_dir <- tempdir()
#' create_evaluation_plots(seurat_obj, results_dir = temp_dir, verbose = FALSE)
#'
#' # Check that plot was created
#' plot_file <- file.path(temp_dir, "umap_experimental.png")
#' file.exists(plot_file)
create_evaluation_plots <- function(seurat_object, results_dir = "./results",
                                    color_palette = NULL, verbose = TRUE) {
  # Create directory if it doesn't exist
  if (!dir.exists(results_dir)) dir.create(results_dir, recursive = TRUE)

  # UMAP plot of experimental data
  p1 <- ggplot2::ggplot(seurat_object@meta.data, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = sample)) +
    ggplot2::geom_point(size = 0.8) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "UMAP of Experimental Data",
      x = "UMAP_1",
      y = "UMAP_2",
      color = "Sample"
    )

  # Save plot
  ggplot2::ggsave(file.path(results_dir, "umap_experimental.png"), p1, width = 10, height = 8)

  # Add more plots as needed...
}
