#' Get best available data layer from Seurat object
#'
#' @description
#' Helper function to determine the best available data layer for FetchData.
#' Checks for normalized data layers and falls back to counts if needed.
#' Handles both Seurat v4 and v5 layer structures.
#'
#' @param seurat_object A Seurat object
#' @return Character string specifying the best layer to use
#' @keywords internal
#' @examples
#' # Example with mock Seurat object layers
#' # This is a simplified example
#' layers <- list(counts = matrix(1:9, nrow = 3), 
#'                data = matrix(2:10, nrow = 3))
#' # In practice, use with actual Seurat object:
#' # best_layer <- get_best_data_layer(seurat_obj)
get_best_data_layer <- function(seurat_object) {
  # Get the default assay (usually "RNA")
  default_assay <- Seurat::DefaultAssay(seurat_object)

  # Check if the assay exists
  if (!default_assay %in% names(seurat_object@assays)) {
    stop(sprintf("Default assay '%s' not found in Seurat object", default_assay))
  }

  # Try to access normalized data layer
  tryCatch(
    {
      data_layer <- seurat_object[[default_assay]]$data
      if (!is.null(data_layer) && nrow(data_layer) > 0) {
        return("data")
      }
    },
    warning = function(w) {
      # Suppress layer warnings and continue to next option
    },
    error = function(e) {
      # Continue to next option if data layer access fails
    }
  )

  # Try to access counts layer as fallback
  tryCatch(
    {
      counts_layer <- seurat_object[[default_assay]]$counts
      if (!is.null(counts_layer) && nrow(counts_layer) > 0) {
        return("counts")
      }
    },
    warning = function(w) {
      # Suppress layer warnings and continue to next option
    },
    error = function(e) {
      # Continue to next option if counts layer access fails
    }
  )

  # Final fallback - let Seurat decide (will likely use default layer)
  return(NULL)
}

#' Parse KEGG pathway file
#'
#' @description
#' Parses a KEGG pathway file (.keg) and extracts gene names for each pathway.
#'
#' @param file_path Character string specifying the path to the KEGG pathway file.
#' Default is "sce00001.keg".
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A list where each element is a character vector of gene names for a pathway.
#' @export
#'
#' @examples
#' \donttest{
#' pathways <- parse_kegg_keg(file_path = "sce00001.keg")
#' }
parse_kegg_keg <- function(file_path = "sce00001.keg", verbose = TRUE) {
  lines <- readLines(file_path)
  pathways <- list()
  current_pathway <- NULL

  for (line in lines) {
    if (startsWith(line, "C")) {
      current_pathway <- sub("^C\\s+\\d+\\s+", "", line)
      pathways[[current_pathway]] <- c()
    } else if (startsWith(line, "D")) {
      matches <- regmatches(line, gregexpr("Y[A-Z][LR][0-9]+[WC]", line))
      if (length(matches[[1]]) > 0) {
        gene_names <- matches[[1]]
        pathways[[current_pathway]] <- c(pathways[[current_pathway]], gene_names)
      }
    }
  }

  return(pathways)
}

#' Build transcriptomic fingerprints
#'
#' @description
#' Builds transcriptomic fingerprints by aggregating gene expression by pathway
#' and creating signature profiles for each condition.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param kegg_pathways A list of KEGG pathways and their associated genes.
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#' \itemize{
#'   \item pathway_matrix: Matrix of pathway expression values
#'   \item signature_matrix: Matrix of condition-specific signature profiles
#' }
#' @export
#'
#' @examples
#' \donttest{
#' fingerprints <- build_fingerprints(seurat_object, kegg_pathways)
#' }
build_fingerprints <- function(seurat_object, kegg_pathways, verbose = TRUE) {
  # Aggregate gene expression by pathway
  pathway_names <- names(kegg_pathways)
  if (verbose) {
    message(sprintf("Processing %d pathways for fingerprint building...", length(pathway_names)))
  }

  # Determine best available data layer
  best_layer <- get_best_data_layer(seurat_object)

  pathway_expression <- sapply(pathway_names, function(pathway) {
    pathway_genes <- kegg_pathways[[pathway]]
    valid_genes <- pathway_genes[!is.na(pathway_genes) & pathway_genes %in% rownames(seurat_object)]
    if (length(valid_genes) == 0) {
      return(rep(NA, ncol(seurat_object)))
    }

    # Use layer parameter only if best_layer is not NULL
    if (is.null(best_layer)) {
      pathway_data <- Seurat::FetchData(seurat_object, vars = valid_genes)
    } else {
      pathway_data <- Seurat::FetchData(seurat_object, vars = valid_genes, layer = best_layer)
    }
    return(rowMeans(pathway_data, na.rm = TRUE))
  }, simplify = FALSE)

  # Create matrix safely with proper column names
  pathway_expression_matrix <- do.call(cbind, pathway_expression)
  colnames(pathway_expression_matrix) <- pathway_names

  # Create signature profiles for each condition
  cell_conditions <- seurat_object@meta.data$sample
  unique_conditions <- unique(cell_conditions)
  signature_profiles <- lapply(unique_conditions, function(condition) {
    condition_rows <- which(cell_conditions == condition)
    colMeans(pathway_expression_matrix[condition_rows, ], na.rm = TRUE)
  })
  signature_matrix <- do.call(cbind, signature_profiles)
  colnames(signature_matrix) <- unique_conditions

  return(list(pathway_matrix = pathway_expression_matrix, signature_matrix = signature_matrix))
}

#' Calculate pathway activities only (for PREDICT mode)
#'
#' @description
#' Calculates pathway activity scores for new data without building signature profiles.
#' This function is used in PREDICT mode where signature matrices are already available
#' from pre-built fingerprints.
#'
#' @details
#' This function performs only the pathway activity calculation part of fingerprint
#' building. It aggregates gene expression by pathway but does not create signature
#' profiles for different conditions. This is the appropriate function to use when
#' applying pre-built fingerprints to new unlabeled data.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param kegg_pathways A list of KEGG pathways with gene names.
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A matrix of pathway activity scores with cells as rows and pathways as columns.
#' @export
#'
#' @examples
#' \donttest{
#' pathway_activities <- calculate_pathway_activities(seurat_object, kegg_pathways)
#' }
calculate_pathway_activities <- function(seurat_object, kegg_pathways, verbose = TRUE) {
  # Aggregate gene expression by pathway
  pathway_names <- names(kegg_pathways)
  if (verbose) {
    message(sprintf("Calculating pathway activities for %d pathways...", length(pathway_names)))
  }

  # Determine best available data layer
  best_layer <- get_best_data_layer(seurat_object)

  pathway_expression <- sapply(pathway_names, function(pathway) {
    pathway_genes <- kegg_pathways[[pathway]]
    valid_genes <- pathway_genes[!is.na(pathway_genes) & pathway_genes %in% rownames(seurat_object)]
    if (length(valid_genes) == 0) {
      return(rep(NA, ncol(seurat_object)))
    }

    # Use layer parameter only if best_layer is not NULL
    if (is.null(best_layer)) {
      pathway_data <- Seurat::FetchData(seurat_object, vars = valid_genes)
    } else {
      pathway_data <- Seurat::FetchData(seurat_object, vars = valid_genes, layer = best_layer)
    }
    return(rowMeans(pathway_data, na.rm = TRUE))
  }, simplify = FALSE)

  # Create matrix safely with proper column names
  pathway_expression_matrix <- do.call(cbind, pathway_expression)
  colnames(pathway_expression_matrix) <- pathway_names

  return(pathway_expression_matrix)
}
