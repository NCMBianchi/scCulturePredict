#' Create UMAP visualization
#'
#' @description
#' Creates a UMAP plot of the single-cell data with customizable aesthetics.
#'
#' @details
#' This function generates a UMAP (Uniform Manifold Approximation and Projection)
#' visualization from a Seurat object. It requires that UMAP coordinates are
#' already computed and stored in the metadata as "UMAP_1" and "UMAP_2".
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters to ensure they meet requirements
#'   \item Checks for the presence of UMAP coordinates in the metadata
#'   \item Creates a ggplot object with UMAP coordinates
#'   \item Customizes the plot appearance based on provided parameters
#'   \item Returns the plot for further customization or direct display
#' }
#'
#' The returned ggplot object can be further customized using standard ggplot2
#' functions and themes. For example, you can add additional layers, change
#' color scales, or modify themes.
#'
#' Error handling ensures that invalid inputs are caught early with informative
#' error messages, making the function robust for interactive and programmatic use.
#'
#' @param seurat_object A Seurat object containing UMAP coordinates in metadata.
#' @param color_by Character string specifying the metadata column to color points by.
#'   Default is "sample".
#' @param point_size Numeric value for point size. Default is 0.8.
#' @param title Character string for plot title. Default is "UMAP Visualization".
#' @param legend_title Character string for legend title. Default is NULL (uses color_by).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A ggplot object containing the UMAP visualization.
#' @export
#'
#' @examples
#' \donttest{
#' plot <- create_umap_plot(seurat_object, color_by = "sample")
#' }
create_umap_plot <- function(seurat_object, color_by = "sample", point_size = 0.8,
                             title = "UMAP Visualization", legend_title = NULL,
                             verbose = TRUE) {
  # Input validation
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object")
  }

  if (!is.character(color_by) || length(color_by) != 1) {
    stop("color_by must be a single character string")
  }

  if (!color_by %in% colnames(seurat_object@meta.data)) {
    stop(sprintf("Column '%s' not found in metadata", color_by))
  }

  if (!is.numeric(point_size) || length(point_size) != 1 || point_size <= 0) {
    stop("point_size must be a single positive number")
  }

  if (!is.null(title) && (!is.character(title) || length(title) != 1)) {
    stop("title must be NULL or a single character string")
  }

  if (!is.null(legend_title) && (!is.character(legend_title) || length(legend_title) != 1)) {
    stop("legend_title must be NULL or a single character string")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Check for UMAP coordinates
  if (!all(c("UMAP_1", "UMAP_2") %in% colnames(seurat_object@meta.data))) {
    stop("UMAP coordinates not found in metadata")
  }

  if (verbose) message("Creating UMAP plot...")

  tryCatch(
    {
      p <- ggplot2::ggplot(seurat_object@meta.data, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = .data[[color_by]])) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = title,
          x = "UMAP_1",
          y = "UMAP_2",
          color = legend_title %||% color_by
        )
      return(p)
    },
    error = function(e) {
      stop(sprintf("Error creating UMAP plot: %s", e$message))
    }
  )
}

#' Create prediction accuracy plot
#'
#' @description
#' Creates a bar plot showing prediction accuracy for different methods and conditions.
#'
#' @details
#' This function visualizes the prediction accuracy results from the
#' \code{\link{evaluate_predictions}} function. It creates a bar plot showing
#' the accuracy of predictions across different samples or conditions.
#'
#' The function processes data differently based on the method parameter:
#' \itemize{
#'   \item "direct": Uses accuracy data from direct similarity-based predictions
#'   \item "svm": Uses accuracy data from Support Vector Machine predictions
#' }
#'
#' The workflow includes:
#' \enumerate{
#'   \item Validating input parameters and checking for required data
#'   \item Extracting the appropriate accuracy data based on the method
#'   \item Creating a bar plot with accuracy percentages
#'   \item Adding percentage labels above each bar
#'   \item Applying custom colors if provided
#' }
#'
#' The resulting plot displays samples on the x-axis and accuracy values on the y-axis,
#' with each bar representing the percentage of correctly predicted cells for that sample.
#' Text labels show the exact percentage values above each bar.
#'
#' The function includes comprehensive error handling to catch missing data,
#' invalid parameters, or processing issues with informative error messages.
#'
#' @param evaluation_results List containing evaluation results from evaluate_predictions().
#' @param method Character string specifying which method to plot ("direct" or "svm").
#' @param title Character string for plot title. Default is NULL.
#' @param color_palette Character vector. Colors to use for the plot
#'   (default: NULL, uses default palette).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A ggplot object containing the accuracy plot.
#' @export
#'
#' @examples
#' \donttest{
#' plot <- create_accuracy_plot(evaluation_results, method = "direct")
#' }
create_accuracy_plot <- function(evaluation_results, method = c("direct", "svm"),
                                 title = NULL, color_palette = NULL, verbose = TRUE) {
  # Input validation
  if (!is.list(evaluation_results)) {
    stop("evaluation_results must be a list")
  }

  method <- match.arg(method)

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

  # Check for required data
  accuracy_data <- if (method == "direct") {
    if (!"direct_accuracy" %in% names(evaluation_results)) {
      stop("direct_accuracy not found in evaluation_results")
    }
    evaluation_results$direct_accuracy
  } else {
    if (!"svm_accuracy" %in% names(evaluation_results)) {
      stop("svm_accuracy not found in evaluation_results")
    }
    evaluation_results$svm_accuracy
  }

  if (verbose) message(sprintf("Creating %s accuracy plot...", method))

  tryCatch(
    {
      p <- ggplot2::ggplot(accuracy_data, ggplot2::aes(x = sample, y = correct, fill = sample)) +
        ggplot2::geom_bar(stat = "identity") +
        ggplot2::geom_text(ggplot2::aes(label = percent), vjust = -0.5) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = title %||% paste(toupper(method), "Prediction Accuracy"),
          x = "Sample",
          y = "Accuracy",
          fill = "Sample"
        ) +
        ggplot2::theme(
          legend.position = "none",
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
        )

      if (!is.null(color_palette)) {
        p <- p + ggplot2::scale_fill_manual(values = color_palette)
      }

      return(p)
    },
    error = function(e) {
      stop(sprintf("Error creating accuracy plot: %s", e$message))
    }
  )
}

#' Create confusion matrix heatmap
#'
#' @description
#' Creates a heatmap visualization of the confusion matrix for prediction results.
#'
#' @details
#' This function creates a heatmap visualization of the confusion matrix from
#' prediction results. A confusion matrix shows the counts of true vs. predicted
#' labels, helping to identify which categories are commonly confused with each other.
#'
#' The function supports three different prediction methods:
#' \itemize{
#'   \item "direct": Uses the confusion matrix from direct similarity-based predictions
#'   \item "threshold": Uses the confusion matrix from threshold-based predictions
#'   \item "svm": Uses the confusion matrix from Support Vector Machine predictions
#' }
#'
#' The workflow consists of:
#' \enumerate{
#'   \item Validating input parameters and checking for required data
#'   \item Extracting the appropriate confusion matrix based on the specified method
#'   \item Converting the matrix to a long-format data frame for plotting
#'   \item Creating a heatmap with tile geometries
#'   \item Adding count labels to each tile
#'   \item Applying a color gradient (customizable via color_palette)
#' }
#'
#' The resulting heatmap displays predicted categories on the x-axis and actual
#' categories on the y-axis. Each cell's color intensity represents the count of
#' predictions, with actual count values displayed as text.
#'
#' The function implements extensive error checking to ensure the required data
#' exists in the evaluation_results object and that all parameters are valid.
#'
#' @param evaluation_results List containing evaluation results from evaluate_predictions().
#' @param method Character string specifying which method to plot ("direct", "threshold", or "svm").
#' @param title Character string for plot title. Default is NULL.
#' @param color_palette Character vector of length 2. Colors for the gradient
#'   (default: NULL, uses white to steelblue).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A ggplot object containing the confusion matrix heatmap.
#' @export
#'
#' @examples
#' \donttest{
#' plot <- create_confusion_heatmap(evaluation_results, method = "svm")
#' }
create_confusion_heatmap <- function(evaluation_results,
                                     method = c("direct", "threshold", "svm"),
                                     title = NULL, color_palette = NULL, verbose = TRUE) {
  # Input validation
  if (!is.list(evaluation_results)) {
    stop("evaluation_results must be a list")
  }

  method <- match.arg(method)

  if (!is.null(title) && (!is.character(title) || length(title) != 1)) {
    stop("title must be NULL or a single character string")
  }

  if (!is.null(color_palette) &&
    (!is.character(color_palette) || length(color_palette) != 2)) {
    stop("color_palette must be NULL or a character vector of length 2")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Check for required data
  confusion_matrix <- switch(method,
    "direct" = {
      if (!"direct_table" %in% names(evaluation_results)) {
        stop("direct_table not found in evaluation_results")
      }
      evaluation_results$direct_table
    },
    "threshold" = {
      if (!"threshold_table" %in% names(evaluation_results)) {
        stop("threshold_table not found in evaluation_results")
      }
      evaluation_results$threshold_table
    },
    "svm" = {
      if (!"svm_table" %in% names(evaluation_results)) {
        stop("svm_table not found in evaluation_results")
      }
      evaluation_results$svm_table
    }
  )

  if (verbose) message(sprintf("Creating %s confusion matrix heatmap...", method))

  tryCatch(
    {
      # Convert to long format for ggplot
      confusion_df <- as.data.frame(confusion_matrix)
      names(confusion_df) <- c("Actual", "Predicted", "Count")

      p <- ggplot2::ggplot(confusion_df, ggplot2::aes(x = Predicted, y = Actual, fill = Count)) +
        ggplot2::geom_tile() +
        ggplot2::geom_text(ggplot2::aes(label = Count), color = "white") +
        ggplot2::scale_fill_gradient(
          low = color_palette[1] %||% "white",
          high = color_palette[2] %||% "steelblue"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::labs(
          title = title %||% paste(toupper(method), "Confusion Matrix"),
          x = "Predicted",
          y = "Actual"
        ) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

      return(p)
    },
    error = function(e) {
      stop(sprintf("Error creating confusion matrix heatmap: %s", e$message))
    }
  )
}

#' Save visualization plots
#'
#' @description
#' Saves multiple visualization plots to a specified directory.
#'
#' @details
#' This function automates the process of creating and saving multiple visualization
#' plots for scCulturePredict analysis results. It generates a complete set of standard
#' visualizations and saves them to the specified directory with consistent naming.
#'
#' The function creates and saves the following plots:
#' \itemize{
#'   \item UMAP plot: Visualization of cells in reduced dimensional space
#'   \item Direct accuracy plot: Bar plot of direct similarity prediction accuracy
#'   \item SVM accuracy plot: Bar plot of SVM prediction accuracy
#'   \item Direct confusion matrix: Heatmap of confusion matrix for direct predictions
#'   \item SVM confusion matrix: Heatmap of confusion matrix for SVM predictions
#' }
#'
#' The workflow consists of:
#' \enumerate{
#'   \item Validating all input parameters
#'   \item Creating the output directory if it doesn't exist
#'   \item Generating each visualization plot using the corresponding functions
#'   \item Saving each plot as a PNG file with the specified prefix
#'   \item Providing progress messages if verbose is TRUE
#' }
#'
#' All plots use consistent formatting and the same color palette (if provided),
#' ensuring visual coherence across the saved visualizations. The prefix parameter
#' allows for organizing multiple analysis runs by giving them distinct filenames.
#'
#' The function includes robust error handling for file system operations and
#' plot creation, with informative error messages if any step fails.
#'
#' @param seurat_object A Seurat object containing analysis results.
#' @param evaluation_results List containing evaluation results from evaluate_predictions().
#' @param output_dir Character string specifying the directory to save plots.
#' @param prefix Character string to prefix plot filenames. Default is "plot".
#' @param color_palette Character vector. Colors to use for the plots
#'   (default: NULL, uses default palette).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return NULL. Plots are saved to the specified directory.
#' @export
#'
#' @examples
#' \donttest{
#' save_visualization_plots(seurat_object, evaluation_results, output_dir = "./plots")
#' }
save_visualization_plots <- function(seurat_object, evaluation_results, output_dir,
                                     prefix = "plot", color_palette = NULL, verbose = TRUE) {
  # Input validation
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object")
  }

  if (!is.list(evaluation_results)) {
    stop("evaluation_results must be a list")
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }

  if (!is.character(prefix) || length(prefix) != 1) {
    stop("prefix must be a single character string")
  }

  if (!is.null(color_palette) &&
    (!is.character(color_palette) || length(color_palette) == 0)) {
    stop("color_palette must be NULL or a non-empty character vector")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Create output directory
  if (verbose) message("Creating output directory...")
  if (!dir.exists(output_dir)) {
    tryCatch(
      {
        dir.create(output_dir, recursive = TRUE)
      },
      error = function(e) {
        stop(sprintf("Error creating output directory: %s", e$message))
      }
    )
  }

  # Create and save plots
  if (verbose) message("Creating and saving plots...")

  tryCatch(
    {
      # UMAP plot
      umap_plot <- create_umap_plot(seurat_object, verbose = verbose)
      ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_umap.png")),
        umap_plot,
        width = 10, height = 8
      )

      # Accuracy plots
      direct_accuracy_plot <- create_accuracy_plot(evaluation_results,
        method = "direct",
        color_palette = color_palette, verbose = verbose
      )
      ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_direct_accuracy.png")),
        direct_accuracy_plot,
        width = 10, height = 8
      )

      svm_accuracy_plot <- create_accuracy_plot(evaluation_results,
        method = "svm",
        color_palette = color_palette, verbose = verbose
      )
      ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_svm_accuracy.png")),
        svm_accuracy_plot,
        width = 10, height = 8
      )

      # Confusion matrix heatmaps
      direct_confusion_plot <- create_confusion_heatmap(evaluation_results,
        method = "direct",
        color_palette = color_palette, verbose = verbose
      )
      ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_direct_confusion.png")),
        direct_confusion_plot,
        width = 10, height = 8
      )

      svm_confusion_plot <- create_confusion_heatmap(evaluation_results,
        method = "svm",
        color_palette = color_palette, verbose = verbose
      )
      ggplot2::ggsave(file.path(output_dir, paste0(prefix, "_svm_confusion.png")),
        svm_confusion_plot,
        width = 10, height = 8
      )

      if (verbose) {
        message("All plots saved successfully")
      }
    },
    error = function(e) {
      stop(sprintf("Error saving plots: %s", e$message))
    }
  )
}
