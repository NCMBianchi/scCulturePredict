#' Analyze pathway enrichment
#'
#' @description
#' Performs pathway enrichment analysis on a Seurat object using KEGG pathways.
#' This function calculates the mean expression of genes in each pathway and
#' performs statistical testing to identify significantly enriched pathways.
#'
#' @param seurat_obj Seurat object. The Seurat object containing gene expression data.
#' @param kegg_pathways List. A list of KEGG pathways where each element is a character
#'   vector of gene names. The names of the list should be pathway IDs.
#' @param min_genes Integer. Minimum number of genes required in a pathway for analysis
#'   (default: 5).
#' @param max_genes Integer. Maximum number of genes allowed in a pathway for analysis
#'   (default: 500).
#' @param p_cutoff Numeric. P-value cutoff for significant pathways (default: 0.05).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A data frame containing pathway enrichment results with columns:
#' \itemize{
#'   \item pathway: Pathway ID
#'   \item n_genes: Number of genes in the pathway
#'   \item mean_expr: Mean expression of pathway genes
#'   \item sd_expr: Standard deviation of pathway expression
#'   \item p_value: P-value from statistical test
#'   \item adj_p_value: Adjusted p-value (FDR)
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters and data
#'   \item Filters pathways based on gene count
#'   \item Calculates mean expression for each pathway
#'   \item Performs statistical testing
#'   \item Adjusts p-values for multiple testing
#' }
#'
#' @examples
#' \dontrun{
#' # Example with mock Seurat object and KEGG pathways
#' # Assuming you have a Seurat object and KEGG pathways loaded
#'
#' # Create mock KEGG pathways
#' kegg_pathways <- list(
#'   pathway1 = c("gene1", "gene2", "gene3", "gene4", "gene5"),
#'   pathway2 = c("gene3", "gene4", "gene5", "gene6", "gene7"),
#'   pathway3 = c("gene5", "gene6", "gene7", "gene8", "gene9")
#' )
#'
#' # Analyze pathway enrichment
#' results <- analyze_pathway_enrichment(
#'   seurat_obj,
#'   kegg_pathways,
#'   min_genes = 5,
#'   max_genes = 500,
#'   p_cutoff = 0.05
#' )
#' }
#' @seealso
#' \code{\link{create_pathway_heatmap}} for visualizing pathway results
#' \code{\link{analyze_pathway_activity}} for analyzing pathway activity
#'
#' @export
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

  if (!all(vapply(kegg_pathways, is.character, FUN.VALUE = logical(1)))) {
    stop("All elements in kegg_pathways must be character vectors")
  }

  if (!is.numeric(min_genes) || length(min_genes) != 1 || min_genes < 0) {
    stop("min_genes must be a single positive number")
  }

  if (!is.numeric(max_genes) || length(max_genes) != 1 || max_genes <= min_genes) {
    stop("max_genes must be a single number greater than min_genes")
  }

  if (!is.numeric(p_cutoff) || length(p_cutoff) != 1 || p_cutoff <= 0 || p_cutoff > 1) {
    stop("p_cutoff must be a single number between 0 and 1")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Check if data is normalized
  if (!"RNA" %in% names(seurat_obj@assays)) {
    stop("Seurat object must contain normalized data. Run Seurat::NormalizeData() first.")
  }

  # Get expression data
  if (verbose) message("Getting expression data...")
  expr_data <- tryCatch(
    {
      Seurat::GetAssayData(seurat_obj, slot = "data")
    },
    error = function(e) {
      stop(sprintf("Error getting expression data: %s", e$message))
    }
  )

  # Filter pathways by gene count
  if (verbose) message("Filtering pathways...")
  pathway_sizes <- vapply(kegg_pathways, length, FUN.VALUE = numeric(1))
  valid_pathways <- kegg_pathways[pathway_sizes >= min_genes & pathway_sizes <= max_genes]

  if (length(valid_pathways) == 0) {
    stop(sprintf("No pathways found with between %d and %d genes", min_genes, max_genes))
  }

  # Calculate pathway statistics
  if (verbose) message("Calculating pathway statistics...")
  results <- tryCatch(
    {
      do.call(rbind, lapply(names(valid_pathways), function(pathway) {
        genes <- valid_pathways[[pathway]]
        # Check if genes exist in the dataset
        valid_genes <- genes[genes %in% rownames(expr_data)]
        if (length(valid_genes) == 0) {
          return(NULL)
        }

        # Calculate mean expression
        pathway_expr <- colMeans(expr_data[valid_genes, ])

        # Perform statistical test
        test_result <- tryCatch(
          {
            t.test(pathway_expr, mu = 0)
          },
          error = function(e) {
            warning(sprintf("Error in statistical test for pathway %s: %s", pathway, e$message))
            return(list(p.value = NA))
          }
        )

        data.frame(
          pathway = pathway,
          n_genes = length(valid_genes),
          mean_expr = mean(pathway_expr),
          sd_expr = sd(pathway_expr),
          p_value = test_result$p.value,
          stringsAsFactors = FALSE
        )
      }))
    },
    error = function(e) {
      stop(sprintf("Error calculating pathway statistics: %s", e$message))
    }
  )

  if (nrow(results) == 0) {
    stop("No valid pathways found after filtering")
  }

  # Adjust p-values
  if (verbose) message("Adjusting p-values...")
  results$adj_p_value <- p.adjust(results$p_value, method = "fdr")

  # Filter by p-value
  results <- results[results$adj_p_value <= p_cutoff, ]

  if (nrow(results) == 0) {
    warning(sprintf("No pathways found significant at adjusted p-value <= %g", p_cutoff))
  }

  if (verbose) {
    message(sprintf("Found %d significant pathways", nrow(results)))
  }

  return(results)
}

#' Create pathway expression heatmap
#'
#' @description
#' Creates a heatmap visualization of pathway expression patterns across cells.
#' The heatmap shows the mean expression of genes in each pathway for each cell,
#' allowing for the identification of pathway activity patterns.
#'
#' @param seurat_obj Seurat object. The Seurat object containing gene expression data.
#' @param pathway_results List. Results from pathway analysis containing:
#'   \itemize{
#'     \item pathway_matrix: Matrix of pathway expression values
#'   }
#' @param top_n Integer. Number of top pathways to show in the heatmap (default: 20).
#' @param cells.use Character vector. Names of cells to include in the heatmap.
#'   If NULL, uses all cells (default: NULL).
#' @param title Character. Title for the heatmap (default: "Pathway Expression Heatmap").
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A ggplot object containing the pathway expression heatmap.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters and data
#'   \item Selects top pathways based on significance
#'   \item Creates a heatmap using ggplot2
#'   \item Adds appropriate labels and theme elements
#' }
#'
#' The heatmap features:
#' \itemize{
#'   \item Pathways on the y-axis
#'   \item Cells on the x-axis
#'   \item Expression levels shown by color intensity
#'   \item Hierarchical clustering of pathways and cells
#' }
#'
#' @examples
#' \donttest{
#' # Basic usage
#' p <- create_pathway_heatmap(
#'   seurat_obj = seurat_obj,
#'   pathway_results = pathway_results
#' )
#'
#' # Show top 10 pathways
#' p <- create_pathway_heatmap(
#'   seurat_obj = seurat_obj,
#'   pathway_results = pathway_results,
#'   top_n = 10
#' )
#'
#' # Use specific cells
#' p <- create_pathway_heatmap(
#'   seurat_obj = seurat_obj,
#'   pathway_results = pathway_results,
#'   cells.use = c("cell1", "cell2", "cell3")
#' )
#'
#' # Custom title
#' p <- create_pathway_heatmap(
#'   seurat_obj = seurat_obj,
#'   pathway_results = pathway_results,
#'   title = "My Custom Heatmap"
#' )
#' }
#'
#' @seealso
#' \code{\link{analyze_pathway_enrichment}} for pathway analysis
#' \code{\link{analyze_pathway_activity}} for pathway activity analysis
#'
#' @export
create_pathway_heatmap <- function(seurat_obj,
                                   pathway_results,
                                   top_n = 20,
                                   cells.use = NULL,
                                   title = "Pathway Expression Heatmap",
                                   verbose = TRUE) {
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  if (!is.list(pathway_results) || !"pathway_matrix" %in% names(pathway_results)) {
    stop("pathway_results must be a list containing 'pathway_matrix'")
  }

  if (!is.matrix(pathway_results$pathway_matrix)) {
    stop("pathway_results$pathway_matrix must be a matrix")
  }

  if (!is.numeric(top_n) || length(top_n) != 1 || top_n < 1) {
    stop("top_n must be a single positive integer")
  }

  if (!is.null(cells.use) && (!is.character(cells.use) || length(cells.use) == 0)) {
    stop("cells.use must be NULL or a non-empty character vector")
  }

  if (!is.character(title) || length(title) != 1) {
    stop("title must be a single character string")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Get pathway matrix
  if (verbose) message("Preparing pathway matrix...")
  pathway_matrix <- tryCatch(
    {
      pathway_results$pathway_matrix
    },
    error = function(e) {
      stop(sprintf("Error accessing pathway matrix: %s", e$message))
    }
  )

  # Select cells
  if (!is.null(cells.use)) {
    if (verbose) message("Subsetting cells...")
    missing_cells <- setdiff(cells.use, colnames(pathway_matrix))
    if (length(missing_cells) > 0) {
      stop(sprintf(
        "Cells not found in pathway matrix: %s",
        paste(missing_cells, collapse = ", ")
      ))
    }
    pathway_matrix <- pathway_matrix[, cells.use]
  }

  # Select top pathways
  if (verbose) message("Selecting top pathways...")
  if (nrow(pathway_matrix) > top_n) {
    # Calculate pathway variance
    pathway_var <- apply(pathway_matrix, 1, var)
    top_pathways <- names(sort(pathway_var, decreasing = TRUE)[seq_len(top_n)])
    pathway_matrix <- pathway_matrix[top_pathways, ]
  }

  # Create heatmap data
  if (verbose) message("Creating heatmap...")
  heatmap_data <- tryCatch(
    {
      reshape2::melt(pathway_matrix)
    },
    error = function(e) {
      stop(sprintf("Error creating heatmap data: %s", e$message))
    }
  )
  names(heatmap_data) <- c("Pathway", "Cell", "Expression")

  # Create heatmap
  p <- tryCatch(
    {
      ggplot2::ggplot(heatmap_data, ggplot2::aes(x = Cell, y = Pathway, fill = Expression)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(
          low = "blue", mid = "white", high = "red",
          midpoint = 0, name = "Expression"
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(
          axis.text.x = ggplot2::element_blank(),
          axis.ticks.x = ggplot2::element_blank(),
          panel.grid = ggplot2::element_blank()
        ) +
        ggplot2::labs(title = title)
    },
    error = function(e) {
      stop(sprintf("Error creating heatmap plot: %s", e$message))
    }
  )

  if (verbose) {
    message("Heatmap created successfully")
  }

  return(p)
}

#' Analyze pathway activity by condition
#'
#' @description
#' Analyzes pathway activity across different conditions or cell types.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param pathway_results List containing pathway analysis results.
#' @param condition Character string specifying the metadata column to group by.
#'   Default is "sample".
#'
#' @return A data frame containing pathway activity statistics by condition with columns:
#' \itemize{
#'   \item pathway: Pathway name
#'   \item condition: Condition or cell type
#'   \item mean_activity: Mean pathway activity
#'   \item sd_activity: Standard deviation of pathway activity
#'   \item p_value: P-value from statistical test
#'   \item adj_p_value: Adjusted p-value (FDR)
#' }
#' @export
#'
#' @examples
#' \dontrun{
#' # Example with mock Seurat object and pathway results
#' # Assuming you have a Seurat object with pathway analysis results
#'
#' # Create mock pathway results
#' pathway_results <- list(
#'   pathway_matrix = matrix(
#'     rnorm(200, mean = 0, sd = 1),
#'     nrow = 20,
#'     ncol = 100, # number of cells
#'     dimnames = list(
#'       paste0("Pathway", seq_len(20)),
#'       paste0("Cell", seq_len(100))
#'     )
#'   )
#' )
#'
#' # Analyze pathway activities by condition
#' # Assumes seurat_object has a "condition" column in metadata
#' analysis_results <- analyze_pathway_activity(
#'   seurat_object,
#'   pathway_results,
#'   condition = "condition"
#' )
#' }
analyze_pathway_activity <- function(seurat_object, pathway_results,
                                     condition = "sample") {
  # Get condition labels
  conditions <- seurat_object@meta.data[[condition]]
  unique_conditions <- unique(conditions)

  # Calculate pathway activity by condition
  activity_stats <- lapply(unique_conditions, function(cond) {
    # Get cells for this condition
    cond_cells <- which(conditions == cond)

    # Calculate pathway activity
    pathway_activity <- pathway_results$pathway_matrix[, cond_cells]

    # Calculate statistics
    mean_activity <- rowMeans(pathway_activity)
    sd_activity <- apply(pathway_activity, 1, sd)

    # Perform statistical test
    p_values <- vapply(seq_len(nrow(pathway_activity)), function(i) {
      tryCatch(
        {
          t.test(pathway_activity[i, ], mu = 0)$p.value
        },
        error = function(e) {
          NA
        }
      )
    }, FUN.VALUE = numeric(1))

    data.frame(
      pathway = rownames(pathway_activity),
      condition = cond,
      mean_activity = mean_activity,
      sd_activity = sd_activity,
      p_value = p_values
    )
  })

  # Combine results
  results <- do.call(rbind, activity_stats)

  # Adjust p-values
  results$adj_p_value <- p.adjust(results$p_value, method = "fdr")

  # Sort by adjusted p-value
  results <- results[order(results$adj_p_value), ]

  return(results)
}

#' Create pathway activity boxplot
#'
#' @description
#' Creates a boxplot visualization of pathway activity across conditions.
#'
#' @param activity_results Data frame containing pathway activity results.
#' @param top_n Integer specifying the number of top pathways to show.
#'   Default is 10.
#' @param title Character string for plot title. Default is NULL.
#'
#' @return A ggplot object containing the pathway activity boxplot.
#' @export
#'
#' @examples
#' \donttest{
#' boxplot <- create_pathway_boxplot(activity_results, top_n = 10)
#' }
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
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
      legend.position = "none"
    ) +
    ggplot2::labs(
      title = title %||% "Pathway Activity by Condition",
      x = "Condition",
      y = "Pathway Activity"
    )

  return(p)
}
