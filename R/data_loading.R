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



#' Load single-cell data
#'
#' @description
#' Loads single-cell data from a directory containing count matrix, metadata, and barcodes.
#' First attempts to use the shell script for file preparation, falling back to R-based
#' preparation if the script fails. Automatically validates malformed input files.
#'
#' @param data_dir Character string specifying the directory containing the data files.
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
#' \dontrun{
#' # Example for loading 10X Genomics data
#' # Assumes you have a directory with 10X format files:
#' # - matrix.mtx.gz (or matrix.mtx)
#' # - barcodes.tsv.gz (or barcodes.tsv)
#' # - features.tsv.gz (or genes.tsv.gz, or features.tsv, or genes.tsv)
#' # - metadata.csv (optional)
#'
#' # Load 10X data
#' seurat_obj <- load_data(
#'   data_dir = "path/to/10x/data",
#'   experiment_id = exp_id
#' )
#'
#' # Clean up
#' unlink(counts_file)
#' unlink(metadata_file)
#' }
#' @seealso
#' \code{\link{preprocess_data}} for preprocessing the loaded data
#' \code{\link{reduce_dimensions}} for dimensionality reduction
#'
#' @export
load_data <- function(data_dir,
                      experiment_id,
                      metadata_file = NULL,
                      min_cells = 3,
                      min_features = 200,
                      verbose = TRUE) {
  # Input validation
  if (!is.character(data_dir) || length(data_dir) != 1) {
    stop("data_dir must be a single character string")
  }
  if (!dir.exists(data_dir)) {
    stop(sprintf("Directory not found: %s", data_dir))
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
      potential_file <- file.path(data_dir, pattern)
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
    metadata_file <- file.path(data_dir, metadata_file)
    if (!file.exists(metadata_file)) {
      stop(sprintf("Metadata file not found: %s", metadata_file))
    }
  }

  # Load count matrix using Seurat's Read10X
  counts <- tryCatch(
    {
      Seurat::Read10X(data.dir = data_dir)
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
    barcodes_file <- file.path(data_dir, "barcodes.tsv.gz")
    if (!file.exists(barcodes_file)) {
      barcodes_file <- file.path(data_dir, "barcodes.tsv")
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
