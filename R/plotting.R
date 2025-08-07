#' Plot scCulture Results
#'
#' @description
#' Automatically generates appropriate visualizations for scCulture analysis results.
#' Detects whether results are from BUILD or PREDICT mode and creates corresponding plots.
#'
#' @param scCulture_results List. Results object returned by scCulture() function.
#'   Must contain seurat_object and mode information.
#' @param plot_type Character. Type of plot to generate. For BUILD mode: "accuracy" (default).
#'   For PREDICT mode: "both" (default), "predictions", or "confidence".
#' @param point_size Numeric. Size of points in UMAP plots (default: 1.5).
#' @param point_alpha Numeric. Transparency of points (default: 0.7 for BUILD, 0.8 for PREDICT).
#' @param return_data Logical. If TRUE, returns the data used for plotting instead of plots (default: FALSE).
#'
#' @return
#' For BUILD mode: A ggplot object showing prediction accuracy on UMAP.
#' For PREDICT mode: A list of ggplot objects containing both predictions and confidence plots (default),
#'   or a single ggplot object when plot_type is "predictions" or "confidence".
#' If return_data = TRUE, returns the data frame used for plotting.
#'
#' @details
#' This function automatically detects the analysis mode from the scCulture results and generates
#' appropriate visualizations:
#'
#' **BUILD Mode Plots:**
#' \itemize{
#'   \item **Accuracy Plot**: UMAP colored by prediction accuracy (blue = correct, red = incorrect)
#'   \item Shows how well the trained model performs on the training data
#'   \item Helps identify regions where the model struggles
#' }
#'
#' **PREDICT Mode Plots:**
#' \itemize{
#'   \item **Both Plots (default)**: Returns both predictions and confidence plots as a list
#'   \item **Predictions Plot**: UMAP colored by predicted culture medium (single plot)
#'   \item **Confidence Plot**: UMAP colored by prediction confidence scores (single plot)
#' }
#'
#' The function extracts UMAP coordinates and relevant metadata from the Seurat object
#' within the scCulture results. It handles missing data gracefully and provides
#' informative error messages for common issues.
#'
#' @examples
#' \dontrun{
#' # Example requires actual scCulture results
#' # Run scCulture first to get results object
#' results <- scCulture(
#'   data_dir = "path/to/data",
#'   mode = "BUILD"
#' )
#'
#' # Create visualization
#' plot <- plot_scCulture(results, plot_type = "umap")
#' print(plot)
#' }
#' @seealso
#' \code{\link{scCulture}} for generating the results object
#'
#' @export
#'
#' @importFrom ggplot2 ggplot aes geom_point scale_color_manual scale_color_viridis_c theme_minimal labs guides guide_legend
#' @importFrom dplyr select mutate
plot_scCulture <- function(scCulture_results,
                           plot_type = NULL,
                           point_size = 1.5,
                           point_alpha = NULL,
                           return_data = FALSE) {
  # Input validation
  if (!is.list(scCulture_results)) {
    stop("scCulture_results must be a list object returned by scCulture()")
  }

  if (!"seurat_object" %in% names(scCulture_results)) {
    stop("scCulture_results must contain a 'seurat_object' component")
  }

  if (!inherits(scCulture_results$seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object")
  }

  if (!is.numeric(point_size) || length(point_size) != 1 || point_size <= 0) {
    stop("point_size must be a single positive number")
  }

  if (!is.null(point_alpha) && (!is.numeric(point_alpha) || length(point_alpha) != 1 ||
    point_alpha < 0 || point_alpha > 1)) {
    stop("point_alpha must be a single number between 0 and 1")
  }

  if (!is.logical(return_data) || length(return_data) != 1) {
    stop("return_data must be a single logical value")
  }

  # Extract Seurat object
  seurat_obj <- scCulture_results$seurat_object

  # Check for required UMAP coordinates
  if (!"UMAP_1" %in% colnames(seurat_obj@meta.data) ||
    !"UMAP_2" %in% colnames(seurat_obj@meta.data)) {
    stop("UMAP coordinates (UMAP_1, UMAP_2) not found in Seurat object metadata")
  }

  # Detect mode from the results object
  mode <- detect_scCulture_mode(scCulture_results)

  # Set default plot_type and point_alpha based on mode
  if (is.null(plot_type)) {
    plot_type <- if (mode == "build") "accuracy" else "both"
  }

  if (is.null(point_alpha)) {
    point_alpha <- if (mode == "build") 0.7 else 0.8
  }

  # Generate plots based on mode
  if (mode == "build") {
    return(create_build_plots(seurat_obj, plot_type, point_size, point_alpha, return_data))
  } else if (mode == "predict") {
    return(create_predict_plots(seurat_obj, plot_type, point_size, point_alpha, return_data))
  } else {
    stop("Unable to determine analysis mode from scCulture results")
  }
}

#' Detect scCulture Analysis Mode
#' @keywords internal
#'
#' @return A character string indicating the detected analysis mode: either 'BUILD' (for training with labeled data) or 'PREDICT' (for prediction on unlabeled data).
#' @examples
#' \dontrun{
#' # This function is internal and requires scCulture results
#' # mode <- detect_scCulture_mode(scCulture_results)
#' }
detect_scCulture_mode <- function(scCulture_results) {
  seurat_obj <- scCulture_results$seurat_object

  # Check for BUILD mode indicators
  has_actual_labels <- "sample" %in% colnames(seurat_obj@meta.data)
  has_predictions <- "classification_pred" %in% colnames(seurat_obj@meta.data)

  # Check for PREDICT mode indicators
  has_fingerprint_source <- "fingerprint_source" %in% names(scCulture_results)
  has_confidence <- "prediction_confidence" %in% colnames(seurat_obj@meta.data)

  # Determine mode based on available data
  if (has_actual_labels && has_predictions && !has_fingerprint_source) {
    return("build")
  } else if (has_predictions && (has_fingerprint_source || has_confidence)) {
    return("predict")
  } else if (has_actual_labels && has_predictions) {
    # Default to build if both labels and predictions present
    return("build")
  } else {
    return("unknown")
    #' @examples
    #' \donttest{
    #' # Example requires a Seurat object with predictions
    #' # plots <- create_build_plots(seurat_obj, output_dir = tempdir())
    #' }
  }
}

#'
#' @return A named list of ggplot2 objects containing visualization plots generated during BUILD mode:
#'   \describe{
#'     \item{accuracy}{Bar plot showing prediction accuracy by sample}
#'     \item{confusion}{Confusion matrix heatmap}
#'     \item{umap}{UMAP plot colored by predictions}
#'     \item{pathway}{Pathway activity heatmap}
#'   }
#' Create BUILD Mode Plots
#' @keywords internal
#' Create plots for BUILD mode
#' @keywords internal
create_build_plots <- function(seurat_obj, plot_type, point_size, point_alpha, return_data) {
  # Validate required columns for BUILD mode
  required_cols <- c("UMAP_1", "UMAP_2", "sample", "classification_pred")
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))

  if (length(missing_cols) > 0) {
    stop(sprintf(
      "BUILD mode requires the following metadata columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  # Create data frame for plotting
  umap_data <- data.frame(
    UMAP_1 = seurat_obj@meta.data$UMAP_1,
    UMAP_2 = seurat_obj@meta.data$UMAP_2,
    Actual = seurat_obj@meta.data$sample,
    Predicted = seurat_obj@meta.data$classification_pred,
    stringsAsFactors = FALSE
  )

  # Calculate prediction accuracy
  umap_data$Correct <- umap_data$Actual == umap_data$Predicted

  # Handle NA predictions
  umap_data$Correct[is.na(umap_data$Predicted)] <- FALSE

  if (return_data) {
    return(umap_data)
  }

  if (plot_type == "accuracy") {
    # Create accuracy plot
    ggplot2::ggplot(umap_data, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = Correct)) +
      ggplot2::geom_point(size = point_size, alpha = point_alpha) +
      ggplot2::scale_color_manual(
        values = c("FALSE" = "red", "TRUE" = "blue"),
        labels = c("FALSE" = "Incorrect", "TRUE" = "Correct")
      ) +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        title = "BUILD Mode: Prediction Accuracy on Training Data",
        subtitle = "Blue = Correct Predictions, Red = Incorrect Predictions",
        x = "UMAP 1",
        #' @examples
        #' \donttest{
        #' # Example requires a Seurat object with predictions
        #' # plots <- create_predict_plots(seurat_obj, output_dir = tempdir())
        #' }
        y = "UMAP 2"
      ) +
      ggplot2::guides(color = ggplot2::guide_legend(title = "Prediction\nAccuracy"))
  } else {
    stop("For BUILD mode, only plot_type = 'accuracy' is supported")
  }
}

#'
#' @return A named list of ggplot2 objects containing visualization plots for PREDICT mode:
#'   \describe{
#'     \item{predictions}{UMAP plot colored by predicted cell types}
#'     \item{confidence}{UMAP plot colored by prediction confidence}
#'     \item{distribution}{Bar plot of predicted cell type frequencies}
#'   }
#' Create PREDICT Mode Plots
#' @keywords internal
#' Create plots for PREDICT mode
#' @keywords internal
create_predict_plots <- function(seurat_obj, plot_type, point_size, point_alpha, return_data) {
  # Validate required columns for PREDICT mode
  base_cols <- c("UMAP_1", "UMAP_2", "classification_pred")
  missing_base <- setdiff(base_cols, colnames(seurat_obj@meta.data))

  if (length(missing_base) > 0) {
    stop(sprintf(
      "PREDICT mode requires the following metadata columns: %s",
      paste(missing_base, collapse = ", ")
    ))
  }

  # Create base data frame
  umap_data <- data.frame(
    UMAP_1 = seurat_obj@meta.data$UMAP_1,
    UMAP_2 = seurat_obj@meta.data$UMAP_2,
    Predicted_Culture = seurat_obj@meta.data$classification_pred,
    stringsAsFactors = FALSE
  )

  # Add confidence if available
  if ("prediction_confidence" %in% colnames(seurat_obj@meta.data)) {
    umap_data$Confidence <- seurat_obj@meta.data$prediction_confidence
  } else {
    umap_data$Confidence <- NA
  }

  if (return_data) {
    return(umap_data)
  }

  # Create plots based on plot_type
  if (plot_type == "predictions") {
    create_predictions_plot(umap_data, point_size, point_alpha)
  } else if (plot_type == "confidence") {
    create_confidence_plot(umap_data, point_size, point_alpha)
  } else if (plot_type == "both") {
    list(
      predictions = create_predictions_plot(umap_data, point_size, point_alpha),
      confidence = create_confidence_plot(umap_data, point_size, point_alpha)
    )
  } else {
    stop("For PREDICT mode, plot_type must be one of: 'predictions', 'confidence', 'both'")
  }
}

#'
#' @return A ggplot2 object visualizing prediction results, typically a UMAP or t-SNE plot with cells colored by their predicted classifications.
#' Create Predictions Plot
#' @keywords internal
#' Create predictions plot
#' @keywords internal
create_predictions_plot <- function(umap_data, point_size, point_alpha) {
  ggplot2::ggplot(umap_data, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = Predicted_Culture)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "PREDICT Mode: Culture Medium Predictions",
      subtitle = "Each cell colored by predicted optimal culture medium",
      x = "UMAP 1",
      y = "UMAP 2"
    ) +
    ggplot2::guides(color = ggplot2::guide_legend(title = "Predicted\nCulture"))
}

#'
#' @return A ggplot2 object showing the distribution of prediction confidence scores. The plot includes a histogram or density plot of confidence values.
#' Create Confidence Plot
#' @keywords internal
#' Create confidence plot
#' @keywords internal
create_confidence_plot <- function(umap_data, point_size, point_alpha) {
  # Check if confidence data is available
  if (all(is.na(umap_data$Confidence))) {
    stop("Confidence scores not available in the results. Cannot create confidence plot.")
  }

  ggplot2::ggplot(umap_data, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = Confidence)) +
    ggplot2::geom_point(size = point_size, alpha = point_alpha) +
    ggplot2::scale_color_viridis_c(name = "Prediction\nConfidence") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "PREDICT Mode: Prediction Confidence Scores",
      subtitle = "Higher values indicate more confident predictions",
      x = "UMAP 1",
      y = "UMAP 2"
    )
}
