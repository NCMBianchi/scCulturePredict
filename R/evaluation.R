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
#' evaluation_results <- evaluate_predictions(seurat_object)
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
#' create_evaluation_plots(seurat_object, results_dir = "./results")
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

#' Evaluate prediction performance
#'
#' @description
#' Evaluates the performance of cell type predictions using various metrics.
#' This function calculates accuracy, precision, recall, F1 score, and confusion
#' matrix for the predictions.
#'
#' @param seurat_obj Seurat object. The Seurat object containing true cell type labels.
#' @param predictions Character vector. Predicted cell type labels.
#' @param true_labels_col Character. Name of the metadata column containing true
#'   cell type labels (default: "cell_type").
#' @param metrics Character vector. Metrics to calculate. Options are "accuracy",
#'   "precision", "recall", "f1", "confusion_matrix" (default: all).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A list containing the requested evaluation metrics:
#' \itemize{
#'   \item accuracy: Overall prediction accuracy
#'   \item precision: Precision for each cell type
#'   \item recall: Recall for each cell type
#'   \item f1: F1 score for each cell type
#'   \item confusion_matrix: Confusion matrix of predictions vs true labels
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters and data
#'   \item Extracts true labels from Seurat object
#'   \item Calculates requested metrics
#'   \item Returns evaluation results
#' }
#'
#' @examples
#' \donttest{
#' # Evaluate all metrics
#' results <- evaluate_predictions(
#'   seurat_obj = seurat_obj,
#'   predictions = predictions
#' )
#'
#' # Evaluate specific metrics
#' results <- evaluate_predictions(
#'   seurat_obj = seurat_obj,
#'   predictions = predictions,
#'   metrics = c("accuracy", "confusion_matrix")
#' )
#' }
#'
#' @seealso
#' \code{\link{predict_cell_types}} for making predictions
#' \code{\link{train_cell_type_classifier}} for training the classifier
#'
#' @export
evaluate_cell_type_predictions <- function(seurat_obj,
                                           predictions,
                                           true_labels_col = "cell_type",
                                           metrics = c("accuracy", "precision", "recall", "f1", "confusion_matrix"),
                                           verbose = TRUE) {
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  if (!is.character(predictions) || length(predictions) != ncol(seurat_obj)) {
    stop("predictions must be a character vector with length matching number of cells")
  }

  if (!is.character(true_labels_col) || length(true_labels_col) != 1) {
    stop("true_labels_col must be a single character string")
  }

  if (!true_labels_col %in% colnames(seurat_obj@meta.data)) {
    stop(sprintf("True labels column '%s' not found in metadata", true_labels_col))
  }

  valid_metrics <- c("accuracy", "precision", "recall", "f1", "confusion_matrix")
  if (!all(metrics %in% valid_metrics)) {
    invalid_metrics <- setdiff(metrics, valid_metrics)
    stop(sprintf(
      "Invalid metrics: %s. Valid options are: %s",
      paste(invalid_metrics, collapse = ", "),
      paste(valid_metrics, collapse = ", ")
    ))
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Get true labels
  if (verbose) message("Getting true labels...")
  true_labels <- tryCatch(
    {
      seurat_obj@meta.data[[true_labels_col]]
    },
    error = function(e) {
      stop(sprintf("Error accessing true labels: %s", e$message))
    }
  )

  if (any(is.na(true_labels))) {
    stop("True labels contain missing values")
  }

  # Calculate metrics
  results <- list()

  if ("confusion_matrix" %in% metrics) {
    if (verbose) message("Calculating confusion matrix...")
    results$confusion_matrix <- tryCatch(
      {
        table(Predicted = predictions, True = true_labels)
      },
      error = function(e) {
        stop(sprintf("Error calculating confusion matrix: %s", e$message))
      }
    )
  }

  if ("accuracy" %in% metrics) {
    if (verbose) message("Calculating accuracy...")
    results$accuracy <- tryCatch(
      {
        mean(predictions == true_labels)
      },
      error = function(e) {
        stop(sprintf("Error calculating accuracy: %s", e$message))
      }
    )
  }

  if (any(c("precision", "recall", "f1") %in% metrics)) {
    if (verbose) message("Calculating per-class metrics...")
    unique_labels <- unique(c(predictions, true_labels))

    if ("precision" %in% metrics) {
      results$precision <- tryCatch(
        {
          sapply(unique_labels, function(label) {
            pred_pos <- predictions == label
            if (sum(pred_pos) == 0) {
              return(0)
            }
            sum(predictions[pred_pos] == true_labels[pred_pos]) / sum(pred_pos)
          })
        },
        error = function(e) {
          stop(sprintf("Error calculating precision: %s", e$message))
        }
      )
    }

    if ("recall" %in% metrics) {
      results$recall <- tryCatch(
        {
          sapply(unique_labels, function(label) {
            true_pos <- true_labels == label
            if (sum(true_pos) == 0) {
              return(0)
            }
            sum(predictions[true_pos] == true_labels[true_pos]) / sum(true_pos)
          })
        },
        error = function(e) {
          stop(sprintf("Error calculating recall: %s", e$message))
        }
      )
    }

    if ("f1" %in% metrics) {
      if (verbose) message("Calculating F1 score...")
      results$f1 <- tryCatch(
        {
          if (!"precision" %in% names(results)) {
            precision <- sapply(unique_labels, function(label) {
              pred_pos <- predictions == label
              if (sum(pred_pos) == 0) {
                return(0)
              }
              sum(predictions[pred_pos] == true_labels[pred_pos]) / sum(pred_pos)
            })
          } else {
            precision <- results$precision
          }

          if (!"recall" %in% names(results)) {
            recall <- sapply(unique_labels, function(label) {
              true_pos <- true_labels == label
              if (sum(true_pos) == 0) {
                return(0)
              }
              sum(predictions[true_pos] == true_labels[true_pos]) / sum(true_pos)
            })
          } else {
            recall <- results$recall
          }

          2 * (precision * recall) / (precision + recall)
        },
        error = function(e) {
          stop(sprintf("Error calculating F1 score: %s", e$message))
        }
      )
    }
  }

  if (verbose) {
    message("Evaluation completed successfully")
  }

  return(results)
}

#' Plot evaluation metrics
#'
#' @description
#' Creates visualizations of prediction evaluation metrics.
#' This function generates various plots to help understand the performance
#' of cell type predictions.
#'
#' @param evaluation_results List. Results from evaluate_predictions().
#' @param plot_type Character. Type of plot to create. Options are:
#'   \itemize{
#'     \item "confusion": Confusion matrix heatmap
#'     \item "metrics": Bar plot of precision, recall, and F1 scores
#'     \item "roc": ROC curve (if probabilities are available)
#'     \item "pr": Precision-recall curve (if probabilities are available)
#'   }
#' @param title Character. Title for the plot (default: NULL).
#' @param color_palette Character vector. Colors to use for the plot
#'   (default: NULL, uses default palette).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A ggplot object containing the requested visualization.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters and data
#'   \item Prepares data for plotting
#'   \item Creates the requested visualization
#'   \item Applies custom styling if specified
#' }
#'
#' @examples
#' \donttest{
#' # Plot confusion matrix
#' p <- plot_evaluation_metrics(
#'   evaluation_results = results,
#'   plot_type = "confusion"
#' )
#'
#' # Plot metrics with custom colors
#' p <- plot_evaluation_metrics(
#'   evaluation_results = results,
#'   plot_type = "metrics",
#'   color_palette = c("red", "blue", "green")
#' )
#'
#' # Plot ROC curve with custom title
#' p <- plot_evaluation_metrics(
#'   evaluation_results = results,
#'   plot_type = "roc",
#'   title = "ROC Curve"
#' )
#' }
#'
#' @details
#' This function creates various types of evaluation plots to visualize prediction performance.
#' It supports four different plot types:
#' \itemize{
#'   \item "confusion": A heatmap visualization of the confusion matrix
#'   \item "metrics": A bar plot of precision, recall, and F1 scores
#'   \item "roc": A receiver operating characteristic curve
#'   \item "pr": A precision-recall curve
#' }
#'
#' The function validates all inputs and checks that the required data for the
#' specified plot type is available in the evaluation_results. It uses ggplot2
#' for visualization and applies appropriate styling and formatting based on the
#' plot type and user-provided parameters.
#'
#' For the confusion matrix plot, cell counts are shown as text labels, while the
#' metrics plot displays precision, recall, and F1 scores with customizable colors.
#' ROC and PR curves visualize classifier performance across different thresholds.
#'
#' @seealso
#' \code{\link{evaluate_predictions}} for calculating evaluation metrics
#'
#' @export
create_evaluation_metrics_plot <- function(evaluation_results,
                                           plot_type = "confusion",
                                           title = NULL,
                                           color_palette = NULL,
                                           verbose = TRUE) {
  # Input validation
  if (!is.list(evaluation_results)) {
    stop("evaluation_results must be a list")
  }

  valid_plot_types <- c("confusion", "metrics", "roc", "pr")
  if (!is.character(plot_type) || length(plot_type) != 1 ||
    !plot_type %in% valid_plot_types) {
    stop(sprintf(
      "plot_type must be one of: %s",
      paste(valid_plot_types, collapse = ", ")
    ))
  }

  if (!is.null(title) && (!is.character(title) || length(title) != 1)) {
    stop("title must be NULL or a single character string")
  }

  if (!is.null(color_palette) &&
    (!is.character(color_palette) || length(color_palette) == 0)) {
    stop("color_palette must be NULL or a non-empty character vector")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Check required data for each plot type
  if (plot_type == "confusion" && !"confusion_matrix" %in% names(evaluation_results)) {
    stop("Confusion matrix not found in evaluation_results")
  }

  if (plot_type == "metrics" &&
    !all(c("precision", "recall", "f1") %in% names(evaluation_results))) {
    stop("Precision, recall, and F1 scores not found in evaluation_results")
  }

  if (plot_type %in% c("roc", "pr") && !"probabilities" %in% names(evaluation_results)) {
    stop("Prediction probabilities not found in evaluation_results")
  }

  # Create plot
  if (verbose) message(sprintf("Creating %s plot...", plot_type))

  p <- tryCatch(
    {
      switch(plot_type,
        "confusion" = {
          # Prepare confusion matrix data
          conf_mat <- evaluation_results$confusion_matrix
          conf_data <- reshape2::melt(conf_mat)
          names(conf_data) <- c("Predicted", "True", "Count")

          # Create heatmap
          ggplot2::ggplot(conf_data, ggplot2::aes(x = True, y = Predicted, fill = Count)) +
            ggplot2::geom_tile() +
            ggplot2::scale_fill_gradient(low = "white", high = "blue") +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
            ggplot2::labs(title = title %||% "Confusion Matrix")
        },
        "metrics" = {
          # Prepare metrics data
          metrics_data <- data.frame(
            CellType = names(evaluation_results$precision),
            Precision = evaluation_results$precision,
            Recall = evaluation_results$recall,
            F1 = evaluation_results$f1
          )
          metrics_data <- reshape2::melt(metrics_data, id.vars = "CellType")

          # Create bar plot
          ggplot2::ggplot(metrics_data, ggplot2::aes(x = CellType, y = value, fill = variable)) +
            ggplot2::geom_bar(stat = "identity", position = "dodge") +
            ggplot2::scale_fill_manual(values = color_palette %||%
              c(
                "Precision" = "blue",
                "Recall" = "red",
                "F1" = "green"
              )) +
            ggplot2::theme_minimal() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
            ggplot2::labs(
              title = title %||% "Prediction Metrics by Cell Type",
              y = "Score",
              fill = "Metric"
            )
        },
        "roc" = {
          # Calculate ROC curve
          roc_data <- pROC::roc(
            evaluation_results$true_labels,
            evaluation_results$probabilities
          )

          # Create ROC curve
          ggplot2::ggplot(data.frame(
            FPR = 1 - roc_data$specificities,
            TPR = roc_data$sensitivities
          ), ggplot2::aes(x = FPR, y = TPR)) +
            ggplot2::geom_line() +
            ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
            ggplot2::theme_minimal() +
            ggplot2::labs(
              title = title %||% "ROC Curve",
              x = "False Positive Rate",
              y = "True Positive Rate"
            )
        },
        "pr" = {
          # Calculate PR curve
          pr_data <- PRROC::pr.curve(
            evaluation_results$true_labels,
            evaluation_results$probabilities
          )

          # Create PR curve
          ggplot2::ggplot(data.frame(
            Recall = pr_data$curve[, 1],
            Precision = pr_data$curve[, 2]
          ), ggplot2::aes(x = Recall, y = Precision)) +
            ggplot2::geom_line() +
            ggplot2::theme_minimal() +
            ggplot2::labs(title = title %||% "Precision-Recall Curve")
        }
      )
    },
    error = function(e) {
      stop(sprintf("Error creating %s plot: %s", plot_type, e$message))
    }
  )

  if (verbose) {
    message("Plot created successfully")
  }

  return(p)
}
