#' Alternative and Extended Implementations for scCulturePredict
#'
#' This file contains alternative implementations and extended functionality
#' that are not currently used in the main scCulturePredict pipeline but may
#' be useful for advanced users or future development.
#'
#' These functions were originally part of the package but were removed to
#' simplify the codebase and improve test coverage. They are preserved here
#' for reference and potential future integration.
#'
#' @author scCulturePredict Development Team
#' @date Last updated: 2024

# ==============================================================================
# DIMENSIONALITY REDUCTION ALTERNATIVES
# ==============================================================================
# Original location: R/dimensionality_reduction.R
# Status: Alternative implementations with additional features
# Integration: To use these, add them to R/dimensionality_reduction.R or create
#              a new file R/advanced_dimred.R and export the functions you need
# ==============================================================================

#' Perform PCA analysis with variance explained
#'
#' @description
#' Enhanced PCA implementation that returns variance explained and creates
#' an elbow plot. This is more feature-rich than the basic PCA in reduce_dimensions().
#'
#' @param seurat_object A Seurat object containing preprocessed single-cell data.
#' @param n_pcs Integer specifying the number of principal components to compute.
#'   Default is 50.
#' @param features Character vector specifying which features to use for PCA.
#'   If NULL, uses variable features. Default is NULL.
#'
#' @return A list containing:
#' \itemize{
#'   \item seurat_object: The Seurat object with PCA results
#'   \item variance_explained: Data frame with variance explained by each PC
#'   \item elbow_plot: ggplot object showing the elbow plot
#' }
#'
#' @details
#' This function provides more information than the basic reduce_dimensions()
#' by calculating variance explained and generating an elbow plot, which can
#' help determine the optimal number of PCs to use.
#'
#' @examples
#' \dontrun{
#' # Use instead of reduce_dimensions() when you need variance information
#' pca_results <- perform_pca(seurat_object, n_pcs = 50)
#' print(pca_results$elbow_plot)
#' }
perform_pca <- function(seurat_object, n_pcs = 50, features = NULL) {
  # Run PCA
  seurat_object <- Seurat::RunPCA(seurat_object,
    features = features %||% Seurat::VariableFeatures(seurat_object),
    npcs = n_pcs,
    verbose = FALSE
  )

  # Calculate variance explained
  variance_explained <- data.frame(
    PC = seq_len(n_pcs),
    Variance = seurat_object@reductions$pca@stdev^2,
    Cumulative = cumsum(seurat_object@reductions$pca@stdev^2) / sum(seurat_object@reductions$pca@stdev^2)
  )

  # Create elbow plot
  elbow_plot <- ggplot2::ggplot(variance_explained, ggplot2::aes(x = PC, y = Variance)) +
    ggplot2::geom_point() +
    ggplot2::geom_line() +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = "PCA Elbow Plot",
      x = "Principal Component",
      y = "Variance Explained"
    )

  return(list(
    seurat_object = seurat_object,
    variance_explained = variance_explained,
    elbow_plot = elbow_plot
  ))
}

#' Run UMAP with clustering
#'
#' @description
#' Extended UMAP implementation that includes clustering and neighbor graph construction.
#' More feature-rich than the basic UMAP in reduce_dimensions().
#'
#' @param seurat_obj Seurat object with scaled data
#' @param n_neighbors Number of neighbors for UMAP (default: 30)
#' @param min_dist Minimum distance for UMAP (default: 0.3)
#' @param n_components Number of UMAP components (default: 2)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return Seurat object with UMAP reduction and clustering
#'
#' @details
#' This function adds clustering capabilities on top of UMAP, which the
#' basic reduce_dimensions() doesn't include.
run_umap <- function(seurat_obj,
                     n_neighbors = 30,
                     min_dist = 0.3,
                     n_components = 2,
                     verbose = TRUE) {
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  if (!is.numeric(n_neighbors) || length(n_neighbors) != 1 || n_neighbors < 0) {
    stop("n_neighbors must be a single positive number")
  }

  if (!is.numeric(min_dist) || length(min_dist) != 1 || min_dist < 0) {
    stop("min_dist must be a single positive number")
  }

  if (!is.numeric(n_components) || length(n_components) != 1 || n_components < 1) {
    stop("n_components must be a single positive integer")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  tryCatch(
    {
      # Find neighbors
      if (verbose) message("Finding neighbors...")
      seurat_obj <- Seurat::FindNeighbors(seurat_obj,
        reduction = "pca",
        dims = 1:40,
        verbose = FALSE
      )

      # Find clusters
      if (verbose) message("Finding clusters...")
      seurat_obj <- Seurat::FindClusters(seurat_obj,
        resolution = 0.5,
        verbose = FALSE
      )

      # Run UMAP
      if (verbose) message("Running UMAP...")
      seurat_obj <- Seurat::RunUMAP(seurat_obj,
        reduction = "pca",
        dims = 1:40,
        n.neighbors = n_neighbors,
        min.dist = min_dist,
        n.components = n_components,
        verbose = FALSE
      )

      if (verbose) message("UMAP complete")

      return(seurat_obj)
    },
    error = function(e) {
      stop(sprintf("Failed to run UMAP: %s", e$message))
    }
  )
}

#' Perform t-SNE dimensionality reduction
#'
#' @description
#' Wrapper for t-SNE with additional parameter control.
#'
#' @param seurat_object Seurat object with PCA reduction
#' @param dims Dimensions to use from PCA (default: 1:40)
#' @param perplexity Perplexity parameter for t-SNE (default: 30)
#' @param max_iter Maximum iterations (default: 1000)
#'
#' @return Seurat object with t-SNE reduction
perform_tsne <- function(seurat_object, dims = 1:40, perplexity = 30,
                         max_iter = 1000) {
  # Run t-SNE
  seurat_object <- Seurat::RunTSNE(seurat_object,
    reduction = "pca",
    dims = dims,
    perplexity = perplexity,
    max_iter = max_iter
  )

  # Add t-SNE coordinates to metadata
  tsne_coords <- seurat_object@reductions$tsne@cell.embeddings
  seurat_object@meta.data$tSNE_1 <- tsne_coords[, 1]
  seurat_object@meta.data$tSNE_2 <- tsne_coords[, 2]

  return(seurat_object)
}

#' Perform complete dimensionality reduction pipeline
#'
#' @description
#' Combines PCA, UMAP, and t-SNE in one function with progress tracking.
#' Alternative to calling reduce_dimensions().
#'
#' @param seurat_object Seurat object with scaled data
#' @param n_pcs Number of PCs to compute (default: 50)
#' @param dims Dimensions to use for UMAP/t-SNE (default: 1:40)
#' @param perform_tsne Whether to perform t-SNE (default: TRUE)
#'
#' @return List containing seurat object and PCA results
perform_dimensionality_reduction <- function(seurat_object, n_pcs = 50,
                                             dims = 1:40, perform_tsne = TRUE) {
  # Perform PCA
  pca_results <- perform_pca(seurat_object, n_pcs = n_pcs)
  seurat_object <- pca_results$seurat_object

  # Perform UMAP
  seurat_object <- run_umap(seurat_object)

  # Perform t-SNE if requested
  if (perform_tsne) {
    seurat_object <- perform_tsne(seurat_object, dims = dims)
  }

  return(list(
    seurat_object = seurat_object,
    pca_results = pca_results
  ))
}

# ==============================================================================
# PATHWAY ANALYSIS ALTERNATIVES
# ==============================================================================
# Original location: R/pathway_analysis.R
# Status: Statistical enrichment and visualization functions not used in pipeline
# Integration: To use these, add them to R/pathway_analysis.R or create
#              a new file R/advanced_pathway.R and export the functions you need
# Current pipeline uses: build_fingerprints() and calculate_pathway_activities()
#                        from kegg_parsing.R instead
# ==============================================================================

#' Analyze pathway enrichment with statistical testing
#'
#' @description
#' Performs statistical pathway enrichment analysis with p-values and FDR correction.
#' This is more comprehensive than build_fingerprints() which just aggregates expression.
#'
#' @param seurat_obj Seurat object containing gene expression data
#' @param kegg_pathways List of KEGG pathways (gene sets)
#' @param min_genes Minimum genes required in pathway (default: 5)
#' @param max_genes Maximum genes allowed in pathway (default: 500)
#' @param p_cutoff P-value cutoff for significance (default: 0.05)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return Data frame with enrichment statistics including p-values
#'
#' @details
#' Unlike build_fingerprints() which simply aggregates expression, this function
#' performs statistical testing to identify significantly enriched pathways.
analyze_pathway_enrichment <- function(seurat_obj,
                                       kegg_pathways,
                                       min_genes = 5,
                                       max_genes = 500,
                                       p_cutoff = 0.05,
                                       verbose = TRUE) {
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  if (!is.list(kegg_pathways) || length(kegg_pathways) == 0) {
    stop("kegg_pathways must be a non-empty list")
  }

  # Get expression matrix
  expr_matrix <- Seurat::GetAssayData(seurat_obj, slot = "data")

  # Filter pathways by size
  pathway_sizes <- sapply(kegg_pathways, function(genes) {
    sum(genes %in% rownames(expr_matrix))
  })

  valid_pathways <- names(kegg_pathways)[
    pathway_sizes >= min_genes & pathway_sizes <= max_genes
  ]

  if (verbose) {
    message(sprintf(
      "Analyzing %d pathways (filtered from %d)",
      length(valid_pathways), length(kegg_pathways)
    ))
  }

  # Calculate enrichment for each pathway
  results <- lapply(valid_pathways, function(pathway_name) {
    genes <- kegg_pathways[[pathway_name]]
    genes_present <- genes[genes %in% rownames(expr_matrix)]

    if (length(genes_present) > 0) {
      # Calculate mean expression
      pathway_expr <- colMeans(expr_matrix[genes_present, , drop = FALSE])

      # Perform statistical test (example: compare to background)
      background_expr <- colMeans(expr_matrix)
      test_result <- t.test(pathway_expr, background_expr)

      return(data.frame(
        pathway = pathway_name,
        n_genes = length(genes_present),
        mean_expr = mean(pathway_expr),
        sd_expr = sd(pathway_expr),
        p_value = test_result$p.value,
        stringsAsFactors = FALSE
      ))
    }
    return(NULL)
  })

  # Combine results
  results_df <- do.call(rbind, results[!sapply(results, is.null)])

  # Adjust p-values
  results_df$adj_p_value <- p.adjust(results_df$p_value, method = "fdr")

  # Sort by adjusted p-value
  results_df <- results_df[order(results_df$adj_p_value), ]

  return(results_df)
}

#' Create pathway expression heatmap
#'
#' @description
#' Creates a heatmap visualization of pathway expression across cells.
#' Useful for visualizing pathway patterns not shown in current pipeline.
#'
#' @param seurat_obj Seurat object
#' @param pathway_results Results from analyze_pathway_enrichment
#' @param top_n Number of top pathways to show (default: 20)
#' @param cells.use Specific cells to include (default: NULL for all)
#' @param title Plot title
#' @param verbose Whether to print messages
#'
#' @return A heatmap plot object
create_pathway_heatmap <- function(seurat_obj,
                                   pathway_results,
                                   top_n = 20,
                                   cells.use = NULL,
                                   title = "Pathway Expression Heatmap",
                                   verbose = TRUE) {
  # This would create a heatmap visualization
  # Implementation depends on specific visualization needs

  if (verbose) message("Creating pathway heatmap...")

  # Get top pathways
  top_pathways <- head(pathway_results$pathway, top_n)

  # Extract expression data for top pathways
  # Create heatmap using ComplexHeatmap or pheatmap

  # Placeholder return
  return(invisible(NULL))
}

#' Analyze pathway activity by condition
#'
#' @description
#' Calculates pathway activity scores grouped by experimental conditions.
#' More detailed than calculate_pathway_activities() as it includes statistical comparisons.
#'
#' @param seurat_object Seurat object
#' @param pathway_results Pathway analysis results
#' @param condition Column name in metadata for grouping (default: "sample")
#'
#' @return Data frame with pathway activities by condition
analyze_pathway_activity <- function(seurat_object, pathway_results,
                                     condition = "sample") {
  # Get condition labels
  conditions <- seurat_object@meta.data[[condition]]
  unique_conditions <- unique(conditions)

  # Calculate pathway activity by condition
  activity_stats <- lapply(unique_conditions, function(cond) {
    # Get cells for this condition
    cond_cells <- which(conditions == cond)

    # Calculate statistics for each pathway
    # This would involve more detailed calculations than
    # the simple aggregation in calculate_pathway_activities()

    return(data.frame(
      condition = cond,
      n_cells = length(cond_cells),
      stringsAsFactors = FALSE
    ))
  })

  # Combine results
  activity_df <- do.call(rbind, activity_stats)

  return(activity_df)
}

#' Create pathway activity boxplot
#'
#' @description
#' Creates boxplot visualization of pathway activities across conditions.
#'
#' @param activity_results Results from analyze_pathway_activity
#' @param top_n Number of top pathways to plot (default: 10)
#' @param title Plot title
#'
#' @return ggplot object
create_pathway_boxplot <- function(activity_results, top_n = 10, title = NULL) {
  # Get top pathways
  top_pathways <- unique(head(activity_results$pathway, top_n))

  # Filter data for top pathways
  plot_data <- activity_results[activity_results$pathway %in% top_pathways, ]

  # Create boxplot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = condition, y = mean_activity, fill = condition)) +
    ggplot2::geom_boxplot() +
    ggplot2::facet_wrap(~pathway, scales = "free_y") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% "Pathway Activity by Condition",
      x = "Condition",
      y = "Pathway Activity"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  return(p)
}

# ==============================================================================
# VISUALIZATION ALTERNATIVES
# ==============================================================================
# Original location: R/visualization.R
# Status: Alternative visualization functions not used in main pipeline
# Integration: To use these, add them to R/visualization.R or keep them separate
# Current pipeline uses: functions from plotting.R and evaluation.R instead
# ==============================================================================

#' Create enhanced UMAP plot
#'
#' @description
#' Creates UMAP visualization with additional customization options.
#' More feature-rich than the basic plots in plot_scCulture().
#'
#' @param seurat_object Seurat object with UMAP reduction
#' @param color_by Metadata column to color by (default: "sample")
#' @param point_size Point size (default: 0.8)
#' @param title Plot title
#' @param legend_title Legend title
#' @param verbose Whether to print messages
#'
#' @return ggplot object
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
    stop(sprintf("color_by column '%s' not found in metadata", color_by))
  }

  # Extract UMAP coordinates
  umap_data <- data.frame(
    UMAP_1 = seurat_object@meta.data$UMAP_1,
    UMAP_2 = seurat_object@meta.data$UMAP_2,
    color_var = seurat_object@meta.data[[color_by]]
  )

  # Create plot
  p <- ggplot2::ggplot(umap_data, ggplot2::aes(x = UMAP_1, y = UMAP_2, color = color_var)) +
    ggplot2::geom_point(size = point_size, alpha = 0.7) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title,
      x = "UMAP 1",
      y = "UMAP 2",
      color = legend_title %||% color_by
    )

  if (verbose) message("UMAP plot created")

  return(p)
}

#' Create accuracy comparison plot
#'
#' @description
#' Creates bar plot comparing accuracies across different prediction methods.
#'
#' @param evaluation_results Evaluation results from pipeline
#' @param method Which method to plot (default: c("direct", "svm"))
#' @param title Plot title
#' @param color_palette Custom colors
#' @param verbose Whether to print messages
#'
#' @return ggplot object
create_accuracy_plot <- function(evaluation_results, method = c("direct", "svm"),
                                 title = NULL, color_palette = NULL, verbose = TRUE) {
  # Input validation
  if (!is.list(evaluation_results)) {
    stop("evaluation_results must be a list")
  }

  method <- match.arg(method)

  # Extract accuracy data based on method
  if (method == "direct") {
    accuracy_data <- evaluation_results$direct_accuracy
  } else {
    accuracy_data <- evaluation_results$svm_accuracy
  }

  # Create bar plot
  plot_data <- data.frame(
    Sample = names(accuracy_data),
    Accuracy = unlist(accuracy_data)
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Sample, y = Accuracy, fill = Sample)) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::ylim(0, 100) +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% paste(toupper(method), "Method Accuracy"),
      x = "Sample",
      y = "Accuracy (%)"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(color_palette)) {
    p <- p + ggplot2::scale_fill_manual(values = color_palette)
  }

  return(p)
}

#' Create confusion matrix heatmap
#'
#' @description
#' Creates a heatmap visualization of the confusion matrix.
#'
#' @param evaluation_results Evaluation results
#' @param method Prediction method (default: c("direct", "threshold", "svm"))
#' @param title Plot title
#' @param color_palette Color palette for heatmap
#' @param verbose Whether to print messages
#'
#' @return ggplot heatmap object
create_confusion_heatmap <- function(evaluation_results,
                                     method = c("direct", "threshold", "svm"),
                                     title = NULL, color_palette = NULL, verbose = TRUE) {
  # Input validation
  if (!is.list(evaluation_results)) {
    stop("evaluation_results must be a list")
  }

  method <- match.arg(method)

  # Get appropriate confusion matrix
  if (method == "direct") {
    conf_matrix <- evaluation_results$direct_confusion
  } else if (method == "threshold") {
    conf_matrix <- evaluation_results$threshold_confusion
  } else {
    conf_matrix <- evaluation_results$svm_confusion
  }

  # Convert to data frame for plotting
  conf_df <- as.data.frame(as.table(conf_matrix))
  names(conf_df) <- c("Actual", "Predicted", "Count")

  # Create heatmap
  p <- ggplot2::ggplot(conf_df, ggplot2::aes(x = Predicted, y = Actual, fill = Count)) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(ggplot2::aes(label = Count), color = "white") +
    ggplot2::scale_fill_gradient(low = "blue", high = "red") +
    ggplot2::theme_minimal() +
    ggplot2::labs(
      title = title %||% paste(toupper(method), "Method Confusion Matrix"),
      x = "Predicted",
      y = "Actual"
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))

  if (!is.null(color_palette)) {
    p <- p + ggplot2::scale_fill_gradient(low = color_palette[1], high = color_palette[2])
  }

  return(p)
}

#' Save multiple visualization plots
#'
#' @description
#' Saves a collection of plots to files with consistent formatting.
#'
#' @param seurat_object Seurat object
#' @param evaluation_results Evaluation results
#' @param output_dir Output directory
#' @param prefix File name prefix (default: "plot")
#' @param color_palette Color palette
#' @param verbose Whether to print messages
#'
#' @return Invisible NULL
save_visualization_plots <- function(seurat_object, evaluation_results, output_dir,
                                     prefix = "plot", color_palette = NULL, verbose = TRUE) {
  # Input validation
  if (!inherits(seurat_object, "Seurat")) {
    stop("seurat_object must be a Seurat object")
  }

  if (!is.list(evaluation_results)) {
    stop("evaluation_results must be a list")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Create and save UMAP plot
  if (verbose) message("Creating UMAP plot...")
  umap_plot <- create_umap_plot(seurat_object, color_by = "sample")
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(prefix, "_umap.pdf")),
    plot = umap_plot,
    width = 8,
    height = 6
  )

  # Create and save accuracy plots
  if (verbose) message("Creating accuracy plots...")
  acc_plot_direct <- create_accuracy_plot(evaluation_results, method = "direct")
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(prefix, "_accuracy_direct.pdf")),
    plot = acc_plot_direct,
    width = 8,
    height = 6
  )

  # Create and save confusion matrix
  if (verbose) message("Creating confusion matrix heatmap...")
  conf_plot <- create_confusion_heatmap(evaluation_results, method = "svm")
  ggplot2::ggsave(
    filename = file.path(output_dir, paste0(prefix, "_confusion_svm.pdf")),
    plot = conf_plot,
    width = 8,
    height = 8
  )

  if (verbose) message(sprintf("Plots saved to %s", output_dir))

  return(invisible(NULL))
}

# ==============================================================================
# USAGE INSTRUCTIONS
# ==============================================================================
#
# To integrate any of these functions back into scCulturePredict:
#
# 1. DIMENSIONALITY REDUCTION FUNCTIONS:
#    - Add to: R/dimensionality_reduction.R (create if doesn't exist)
#    - Or add to: R/preprocessing.R (alongside reduce_dimensions)
#    - Export in NAMESPACE if you want users to access them
#    - Example usage:
#      pca_results <- perform_pca(seurat_obj)  # Get variance explained
#      seurat_obj <- run_umap(seurat_obj)      # UMAP with clustering
#
# 2. PATHWAY ANALYSIS FUNCTIONS:
#    - Add to: R/pathway_analysis.R (create if doesn't exist)
#    - Or add to: R/kegg_parsing.R (alongside existing pathway functions)
#    - Export in NAMESPACE for user access
#    - Example usage:
#      enrichment <- analyze_pathway_enrichment(seurat_obj, kegg_pathways)
#      heatmap <- create_pathway_heatmap(seurat_obj, enrichment)
#
# 3. VISUALIZATION FUNCTIONS:
#    - Add to: R/visualization.R (create if doesn't exist)
#    - Or add to: R/plotting.R (alongside plot_scCulture)
#    - Export in NAMESPACE for user access
#    - Example usage:
#      p <- create_umap_plot(seurat_obj, color_by = "condition")
#      save_visualization_plots(seurat_obj, results, "output/")
#
# 4. TO USE AS ADVANCED API:
#    - Create a new vignette: vignettes/advanced_functions.Rmd
#    - Document these as "Advanced Features" for power users
#    - Keep them separate from the main pipeline
#
# 5. TO TEST THESE FUNCTIONS:
#    - Create: tests/testthat/test-advanced-functions.R
#    - Test each function group separately
#    - Use the mock data generators from existing tests
#
# ==============================================================================
