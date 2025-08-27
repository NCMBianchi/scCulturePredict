#' Validate and fix malformed input files
#'
#' @description
#' Checks for malformed input files (e.g., files with "x" in first row) and fixes them.
#'
#' @param file_path Character string specifying the path to the file to validate
#' @param sep Character string specifying the field separator
#' @param header Logical indicating whether the file has a header
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A data frame with the corrected data
#' @keywords internal
validate_and_fix_file <- function(file_path, sep = "\t", header = TRUE, verbose = TRUE) {
  # Input validation
  if (!is.character(file_path) || length(file_path) != 1) {
    stop("file_path must be a single character string")
  }
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }
  if (!is.character(sep) || length(sep) != 1) {
    stop("sep must be a single character string")
  }
  if (!is.logical(header) || length(header) != 1) {
    stop("header must be a single logical value")
  }
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  tryCatch(
    {
      # Read the first few lines to check for malformed data
      if (verbose) message("Reading file header...")
      lines <- readLines(file_path, n = 5)

      # Check if first line contains only 'x' or similar malformed data
      if (grepl("^x\\s*$", lines[1]) || grepl("^x\\s+.*$", lines[1])) {
        if (verbose) message("Detected malformed header, skipping first line...")
        # Read the file skipping the first line
        data <- read.csv(file_path, sep = sep, header = header, skip = 1)
        return(data)
      }

      # If file is not malformed, read normally
      if (verbose) message("Reading file...")
      return(read.csv(file_path, sep = sep, header = header))
    },
    error = function(e) {
      stop(sprintf("Failed to read file %s: %s", file_path, e$message))
    }
  )
}



#' Load 10X Genomics single-cell data
#'
#' @description
#' Loads single-cell data from a directory containing count matrix, metadata, and barcodes.
#' First attempts to use the shell script for file preparation, falling back to R-based
#' preparation if the script fails. Automatically validates malformed input files.
#'
#' @param tenx_data_dir Character string specifying the directory containing the data files.
#'   The directory should contain files with the following naming pattern:
#'   \itemize{
#'     \item \code{<experiment_id>_metadata.csv} - Cell metadata file
#'     \item \code{<experiment_id>_counts.csv} - Gene expression count matrix
#'   }
#' @param experiment_id Character string specifying the experiment ID prefix in filenames.
#'   This should match the prefix used in your data files.
#' @param metadata_file Optional. Character string specifying the name of the metadata file.
#'   If NULL, will be automatically detected based on experiment_id.
#' @param min_cells Integer specifying the minimum number of cells expressing a gene
#'   for the gene to be included. Default is 3.
#' @param min_features Integer specifying the minimum number of genes expressed in a cell
#'   for the cell to be included. Default is 200.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#'
#' @return A Seurat object containing:
#' \itemize{
#'   \item Count matrix in the \code{counts} slot
#'   \item Cell metadata in the \code{meta.data} slot
#'   \item Filtered genes and cells based on \code{min_cells} and \code{min_features}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters and file existence
#'   \item Automatically detects metadata and counts files if not specified
#'   \item Loads and validates metadata and counts data
#'   \item Checks for data consistency (matching cell numbers)
#'   \item Creates a Seurat object with appropriate filtering
#' }
#'
#' @examples
#' # Load 10X data using package's built-in example data
#' data_dir <- system.file("extdata", "example_data_10x",
#'   package = "scCulturePredict"
#' )
#'
#' # Check that example data exists
#' if (dir.exists(data_dir) && length(list.files(data_dir)) > 0) {
#'   # Load the 10X format data
#'   seurat_obj <- load_10x_data(
#'     tenx_data_dir = data_dir,
#'     experiment_id = "example_10x",
#'     min_cells = 3,
#'     min_features = 10,
#'     verbose = FALSE
#'   )
#'
#'   # Display basic information about the loaded object
#'   print(seurat_obj)
#'   cat("Number of cells:", ncol(seurat_obj), "\n")
#'   cat("Number of genes:", nrow(seurat_obj), "\n")
#' }
#' @seealso
#' \code{\link{preprocess_data}} for preprocessing the loaded data
#' \code{\link{reduce_dimensions}} for dimensionality reduction
#'
#' @export
load_10x_data <- function(tenx_data_dir,
                          experiment_id,
                          metadata_file = NULL,
                          min_cells = 3,
                          min_features = 200,
                          verbose = TRUE) {
  # Input validation
  if (!is.character(tenx_data_dir) || length(tenx_data_dir) != 1) {
    stop("tenx_data_dir must be a single character string")
  }
  if (!dir.exists(tenx_data_dir)) {
    stop(sprintf("Directory not found: %s", tenx_data_dir))
  }

  if (!is.character(experiment_id) || length(experiment_id) != 1) {
    stop("experiment_id must be a single character string")
  }

  if (!is.null(metadata_file) && (!is.character(metadata_file) || length(metadata_file) != 1)) {
    stop("metadata_file must be a single character string or NULL")
  }

  if (!is.numeric(min_cells) || length(min_cells) != 1 || min_cells < 0) {
    stop("min_cells must be a single positive number")
  }

  if (!is.numeric(min_features) || length(min_features) != 1 || min_features < 0) {
    stop("min_features must be a single positive number")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Load 10X Genomics format data
  if (verbose) message("Loading 10X Genomics data...")

  # Find metadata file for 10X format
  if (is.null(metadata_file)) {
    metadata_patterns <- c(
      paste0(experiment_id, "_metadata.tsv.gz"),
      paste0(experiment_id, "_metadata.tsv"),
      "metadata.tsv.gz",
      "metadata.tsv"
    )

    for (pattern in metadata_patterns) {
      potential_file <- file.path(tenx_data_dir, pattern)
      if (file.exists(potential_file)) {
        metadata_file <- potential_file
        break
      }
    }

    if (is.null(metadata_file)) {
      if (verbose) message("No metadata file found, will create basic metadata from barcodes")
      metadata_file <- NULL
    }
  } else {
    metadata_file <- file.path(tenx_data_dir, metadata_file)
    if (!file.exists(metadata_file)) {
      stop(sprintf("Metadata file not found: %s", metadata_file))
    }
  }

  # Load count matrix using Seurat's Read10X
  counts <- tryCatch(
    {
      Seurat::Read10X(data.dir = tenx_data_dir)
    },
    error = function(e) {
      stop(sprintf("Error loading 10X data: %s", e$message))
    }
  )

  # Create initial Seurat object without metadata first
  seurat_obj <- tryCatch(
    {
      Seurat::CreateSeuratObject(
        counts = counts,
        min.cells = min_cells,
        min.features = min_features
      )
    },
    error = function(e) {
      stop(sprintf("Error creating Seurat object: %s", e$message))
    }
  )

  # Load and add metadata if available (using original approach)
  if (!is.null(metadata_file)) {
    if (verbose) message("Loading metadata...")

    # Load metadata without row names (like original script)
    metadata <- tryCatch(
      {
        if (grepl("\\.csv$", metadata_file)) {
          read.csv(metadata_file, row.names = NULL)
        } else {
          read.table(metadata_file, header = TRUE, sep = "\t", row.names = NULL)
        }
      },
      error = function(e) {
        stop(sprintf("Error loading metadata file: %s", e$message))
      }
    )

    # Rename first column to 'barcode' (like original script)
    colnames(metadata)[1] <- "barcode"

    # Load barcodes file
    barcodes_file <- file.path(tenx_data_dir, "barcodes.tsv.gz")
    if (!file.exists(barcodes_file)) {
      barcodes_file <- file.path(tenx_data_dir, "barcodes.tsv")
    }

    if (file.exists(barcodes_file)) {
      barcodes <- tryCatch(
        {
          read.table(barcodes_file, sep = "\t", header = FALSE, col.names = c("index", "barcode"))
        },
        error = function(e) {
          # If barcodes file doesn't exist or fails, use existing column names
          if (verbose) message("Using existing cell names as barcodes")
          barcodes <- data.frame(index = seq_len(ncol(seurat_obj)), barcode = colnames(seurat_obj))
        }
      )

      # Ensure cell names match barcodes (like original script)
      barcode_match <- match(barcodes$barcode, colnames(seurat_obj))
      if (!any(is.na(barcode_match))) {
        colnames(seurat_obj) <- barcodes$barcode[barcode_match]
      }
    }

    # Validate sample column exists
    if (!"sample" %in% colnames(metadata)) {
      stop("Metadata must contain a 'sample' column")
    }

    # Add metadata using AddMetaData (like original script)
    metadata$cell <- rownames(metadata)
    for (meta_column in setdiff(colnames(metadata), c("barcode", "index", "cell"))) {
      seurat_obj <- Seurat::AddMetaData(seurat_obj,
        metadata = metadata[[meta_column]],
        col.name = meta_column
      )
    }

    # Rename existing columns if needed (like original script)
    existing_cols <- intersect(
      c("PC_1", "PC_2", "PC_3", "UMAP_1", "UMAP_2"),
      colnames(seurat_obj@meta.data)
    )
    for (col_name in existing_cols) {
      new_col_name <- paste0(col_name, "_old")
      names(seurat_obj@meta.data)[names(seurat_obj@meta.data) == col_name] <- new_col_name
    }
  } else {
    # Create basic sample metadata if none provided
    if (verbose) message("Creating basic metadata from cell barcodes...")
    seurat_obj@meta.data$sample <- "sample1"
  }

  if (verbose) {
    message(sprintf(
      "Created Seurat object with %d cells and %d genes",
      ncol(seurat_obj), nrow(seurat_obj)
    ))
  }

  return(seurat_obj)
}

#' Load SingleCellExperiment data
#'
#' @description
#' Loads single-cell data from a SingleCellExperiment object or RDS file containing one.
#' Converts the SingleCellExperiment to a Seurat object for pipeline compatibility.
#'
#' @param sce_data_path A SingleCellExperiment object or character string path to an RDS file
#'   containing a SingleCellExperiment object.
#' @param experiment_id Character string specifying the experiment ID for tracking.
#' @param min_cells Integer specifying the minimum number of cells expressing a gene
#'   for the gene to be included. Default is 3.
#' @param min_features Integer specifying the minimum number of genes expressed in a cell
#'   for the cell to be included. Default is 200.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#' @param handle_duplicates Character string specifying how to handle duplicate gene names.
#'   Options are "make_unique" (default), "aggregate", "first", or "error".
#'   "make_unique" appends .1, .2 etc to duplicates, "aggregate" sums duplicate genes,
#'   "first" keeps only the first occurrence, "error" stops with an informative message.
#'
#' @return A Seurat object containing:
#' \itemize{
#'   \item Count matrix from the SingleCellExperiment
#'   \item Cell metadata from colData
#'   \item Filtered genes and cells based on \code{min_cells} and \code{min_features}
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Loads SCE object from file if path is provided
#'   \item Validates the SingleCellExperiment object structure
#'   \item Extracts counts and metadata
#'   \item Converts to Seurat object with appropriate filtering
#'   \item Preserves sample information from colData
#' }
#'
#' @examples
#' # Create a small example SingleCellExperiment
#' if (requireNamespace("SingleCellExperiment", quietly = TRUE)) {
#'   library(SingleCellExperiment)
#'
#'   # Generate example count data
#'   set.seed(123)
#'   counts <- matrix(rpois(2000, lambda = 5), ncol = 40, nrow = 50)
#'   colnames(counts) <- paste0("Cell", seq_len(40))
#'   rownames(counts) <- paste0("Gene", seq_len(50))
#'
#'   # Create SingleCellExperiment object
#'   sce <- SingleCellExperiment(assays = list(counts = counts))
#'
#'   # Add sample metadata
#'   colData(sce)$sample <- rep(c("SampleA", "SampleB"), each = 20)
#'
#'   # Load the SCE object into Seurat format
#'   seurat_obj <- load_sce_data(
#'     sce_data_path = sce,
#'     experiment_id = "example_sce",
#'     min_cells = 3,
#'     min_features = 10,
#'     verbose = FALSE
#'   )
#'
#'   # Display basic information
#'   print(seurat_obj)
#'   cat("Number of cells:", ncol(seurat_obj), "\n")
#'   cat("Number of genes:", nrow(seurat_obj), "\n")
#' }
#'
#' @seealso
#' \code{\link{load_10x_data}} for loading 10X Genomics format data
#' \code{\link{preprocess_data}} for preprocessing the loaded data
#'
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom SummarizedExperiment assay colData
#' @export
load_sce_data <- function(sce_data_path,
                          experiment_id,
                          min_cells = 3,
                          min_features = 200,
                          verbose = TRUE,
                          handle_duplicates = "make_unique") {
  # Input validation
  if (!is.character(experiment_id) || length(experiment_id) != 1) {
    stop("experiment_id must be a single character string")
  }

  if (!is.numeric(min_cells) || length(min_cells) != 1 || min_cells < 0) {
    stop("min_cells must be a single positive number")
  }

  if (!is.numeric(min_features) || length(min_features) != 1 || min_features < 0) {
    stop("min_features must be a single positive number")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  if (!is.character(handle_duplicates) || length(handle_duplicates) != 1) {
    stop("handle_duplicates must be a single character string")
  }

  valid_duplicate_methods <- c("make_unique", "aggregate", "first", "error")
  if (!handle_duplicates %in% valid_duplicate_methods) {
    stop(sprintf("handle_duplicates must be one of: %s", paste(valid_duplicate_methods, collapse = ", ")))
  }

  # Load SCE from file if path is provided
  if (is.character(sce_data_path)) {
    if (length(sce_data_path) != 1) {
      stop("sce_data_path must be a single character string")
    }
    if (!file.exists(sce_data_path)) {
      stop(sprintf("File not found: %s", sce_data_path))
    }

    if (verbose) message(sprintf("Loading SingleCellExperiment from %s...", sce_data_path))

    sce_data_path <- tryCatch(
      {
        readRDS(sce_data_path)
      },
      error = function(e) {
        stop(sprintf("Failed to load RDS file: %s", e$message))
      }
    )
  }

  # Validate it's a SingleCellExperiment
  if (!inherits(sce_data_path, "SingleCellExperiment")) {
    stop("sce_data_path must be a SingleCellExperiment object or path to RDS file containing one")
  }

  if (verbose) message("Converting SingleCellExperiment to Seurat object...")

  # Extract counts matrix
  # Try different assay names commonly used
  counts <- NULL
  for (assay_name in c("counts", "raw", "rawcounts")) {
    if (assay_name %in% names(SummarizedExperiment::assays(sce_data_path))) {
      counts <- SummarizedExperiment::assay(sce_data_path, assay_name)
      if (verbose) message(sprintf("Using '%s' assay for count data", assay_name))
      break
    }
  }

  if (is.null(counts)) {
    # If no counts found, use the first assay
    if (length(SummarizedExperiment::assays(sce_data_path)) > 0) {
      counts <- SummarizedExperiment::assay(sce_data_path, 1)
      if (verbose) message("Using first assay as count data")
    } else {
      stop("No assays found in SingleCellExperiment object")
    }
  }

  # Ensure counts matrix has proper rownames and colnames
  if (is.null(rownames(counts))) {
    # Get row names from SCE object
    rownames(counts) <- rownames(sce_data_path)
    if (verbose) message("Added row names from SingleCellExperiment object")
  }

  if (is.null(colnames(counts))) {
    # Get column names from SCE object
    colnames(counts) <- colnames(sce_data_path)
    if (verbose) message("Added column names from SingleCellExperiment object")
  }

  # Check for empty matrix before setting names
  if (nrow(counts) == 0 || ncol(counts) == 0) {
    stop("Cannot load SingleCellExperiment with empty count matrix (0 genes or 0 cells)")
  }

  # Final validation
  if (is.null(rownames(counts)) || length(rownames(counts)) == 0) {
    # If still no rownames, create generic ones
    rownames(counts) <- paste0("Feature", seq_len(nrow(counts)))
    if (verbose) message("Created generic feature names")
  }

  if (is.null(colnames(counts)) || length(colnames(counts)) == 0) {
    # If still no colnames, create generic ones
    colnames(counts) <- paste0("Cell", seq_len(ncol(counts)))
    if (verbose) message("Created generic cell names")
  }

  # Handle duplicate gene names before creating Seurat object
  if (any(duplicated(rownames(counts)))) {
    n_dups <- sum(duplicated(rownames(counts)))
    dup_genes <- unique(rownames(counts)[duplicated(rownames(counts))])

    if (handle_duplicates == "error") {
      stop(sprintf(
        "Found %d duplicate gene names: %s. Please resolve before loading or use a different handle_duplicates option.",
        n_dups, paste(head(dup_genes, 10), collapse = ", ")
      ))
    }

    if (verbose) {
      message(sprintf("Found %d duplicate gene names, handling with method: %s", n_dups, handle_duplicates))
    }

    if (handle_duplicates == "make_unique") {
      rownames(counts) <- make.unique(rownames(counts))
      if (verbose) message("Made gene names unique by appending .1, .2, etc. to duplicates")
    } else if (handle_duplicates == "first") {
      counts <- counts[!duplicated(rownames(counts)), ]
      if (verbose) message(sprintf("Kept only first occurrence of duplicate genes, removed %d rows", n_dups))
    } else if (handle_duplicates == "aggregate") {
      # Sum duplicate genes
      if (verbose) message("Aggregating duplicate genes by summing counts...")
      unique_genes <- unique(rownames(counts))
      aggregated_counts <- matrix(0, nrow = length(unique_genes), ncol = ncol(counts))
      rownames(aggregated_counts) <- unique_genes
      colnames(aggregated_counts) <- colnames(counts)

      for (i in seq_along(unique_genes)) {
        gene <- unique_genes[i]
        gene_rows <- which(rownames(counts) == gene)
        if (length(gene_rows) == 1) {
          aggregated_counts[i, ] <- counts[gene_rows, ]
        } else {
          aggregated_counts[i, ] <- colSums(counts[gene_rows, , drop = FALSE])
        }
      }
      counts <- aggregated_counts
      if (verbose) message(sprintf("Aggregated %d duplicate genes", n_dups))
    }
  }

  # Create Seurat object
  seurat_obj <- tryCatch(
    {
      Seurat::CreateSeuratObject(
        counts = counts,
        min.cells = min_cells,
        min.features = min_features,
        project = experiment_id
      )
    },
    error = function(e) {
      stop(sprintf("Error creating Seurat object: %s", e$message))
    }
  )

  # Extract and add metadata from colData
  col_data <- as.data.frame(SummarizedExperiment::colData(sce_data_path))

  if (ncol(col_data) > 0) {
    if (verbose) message("Adding metadata from colData...")

    # Get cells that passed filtering
    filtered_cells <- colnames(seurat_obj)

    # Subset col_data to only include cells that are in the Seurat object
    # This handles the case where cells were filtered out during CreateSeuratObject
    if (all(filtered_cells %in% rownames(col_data))) {
      # Cell names match, subset by name
      col_data <- col_data[filtered_cells, , drop = FALSE]
      if (verbose) message(sprintf("Subset metadata to %d cells that passed filtering", length(filtered_cells)))
    } else if (all(filtered_cells %in% colnames(counts))) {
      # Try to match using original column names from counts matrix
      original_cells <- colnames(counts)
      kept_indices <- which(original_cells %in% filtered_cells)
      col_data <- col_data[kept_indices, , drop = FALSE]
      rownames(col_data) <- filtered_cells
      if (verbose) message(sprintf("Matched metadata by position for %d filtered cells", length(filtered_cells)))
    } else {
      # Last resort: if the number of cells matches exactly, assume they're in order
      if (ncol(seurat_obj) == nrow(col_data)) {
        rownames(col_data) <- filtered_cells
        if (verbose) message("Aligned metadata rows with Seurat cell names by position")
      } else {
        # Only some cells were kept, but we can't match them
        if (verbose) {
          message(sprintf(
            "Warning: Cannot match metadata. Seurat has %d cells, metadata has %d rows",
            ncol(seurat_obj), nrow(col_data)
          ))
          message("Skipping metadata addition")
        }
        col_data <- data.frame(row.names = filtered_cells)
      }
    }

    # Add each metadata column
    for (col_name in colnames(col_data)) {
      # Skip if column already exists in Seurat metadata
      if (col_name %in% colnames(seurat_obj@meta.data)) {
        if (verbose) message(sprintf("Skipping existing metadata column: %s", col_name))
        next
      }

      seurat_obj <- Seurat::AddMetaData(
        seurat_obj,
        metadata = col_data[[col_name]],
        col.name = col_name
      )
    }

    # Ensure 'sample' column exists
    if (!"sample" %in% colnames(seurat_obj@meta.data)) {
      if ("Sample" %in% colnames(seurat_obj@meta.data)) {
        # Rename 'Sample' to 'sample' if it exists
        seurat_obj@meta.data$sample <- seurat_obj@meta.data$Sample
        if (verbose) message("Renamed 'Sample' column to 'sample'")
      } else {
        # Create default sample column
        if (verbose) message("No 'sample' column found, creating default...")
        seurat_obj@meta.data$sample <- "sample1"
      }
    }
  } else {
    # No metadata in colData, create basic metadata
    if (verbose) message("No metadata in colData, creating basic metadata...")
    seurat_obj@meta.data$sample <- "sample1"
  }

  # Add reduced dimensions if they exist
  if (length(SingleCellExperiment::reducedDims(sce_data_path)) > 0) {
    if (verbose) message("Transferring dimensionality reductions...")

    red_dims <- SingleCellExperiment::reducedDims(sce_data_path)

    # Transfer PCA if exists
    if ("PCA" %in% names(red_dims)) {
      pca_data <- red_dims[["PCA"]]
      # Note: Would need additional code to properly transfer as Seurat DimReduc object
      # For now, add to metadata as coordinates
      if (ncol(pca_data) >= 2) {
        seurat_obj@meta.data$PC_1_sce <- pca_data[, 1]
        seurat_obj@meta.data$PC_2_sce <- pca_data[, 2]
      }
    }

    # Transfer UMAP if exists
    if ("UMAP" %in% names(red_dims)) {
      umap_data <- red_dims[["UMAP"]]
      if (ncol(umap_data) >= 2) {
        seurat_obj@meta.data$UMAP_1_sce <- umap_data[, 1]
        seurat_obj@meta.data$UMAP_2_sce <- umap_data[, 2]
      }
    }

    # Transfer TSNE if exists
    if ("TSNE" %in% names(red_dims)) {
      tsne_data <- red_dims[["TSNE"]]
      if (ncol(tsne_data) >= 2) {
        seurat_obj@meta.data$tSNE_1_sce <- tsne_data[, 1]
        seurat_obj@meta.data$tSNE_2_sce <- tsne_data[, 2]
      }
    }
  }

  if (verbose) {
    message(sprintf(
      "Created Seurat object with %d cells and %d genes from SingleCellExperiment",
      ncol(seurat_obj), nrow(seurat_obj)
    ))
  }

  return(seurat_obj)
}

#' Process metadata for Seurat object
#'
#' @description
#' Internal helper function to process metadata and barcodes for a Seurat object
#'
#' @param seurat_object A Seurat object
#' @param metadata Data frame containing metadata
#' @param barcodes Data frame containing barcodes
#'
#' @return A Seurat object with processed metadata
#' @keywords internal
process_metadata <- function(seurat_object, metadata, barcodes) {
  # Ensure cell names match barcodes
  barcode_match <- match(barcodes$barcode, colnames(seurat_object))
  if (any(is.na(barcode_match))) {
    warning("Some barcodes do not match cells in the Seurat object.")
  } else {
    colnames(seurat_object) <- barcodes$barcode[barcode_match]
  }

  # Add metadata to Seurat object
  metadata$cell <- rownames(metadata)
  for (meta_column in setdiff(colnames(metadata), c("barcode", "index", "cell"))) {
    seurat_object <- AddMetaData(seurat_object, metadata = metadata[[meta_column]], col.name = meta_column)
  }

  # Rename existing columns if needed
  existing_cols <- intersect(c("PC_1", "PC_2", "PC_3", "UMAP_1", "UMAP_2"), colnames(seurat_object@meta.data))
  for (col_name in existing_cols) {
    new_col_name <- paste0(col_name, "_old")
    names(seurat_object@meta.data)[names(seurat_object@meta.data) == col_name] <- new_col_name
  }

  return(seurat_object)
}
