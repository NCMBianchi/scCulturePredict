#' Load required packages for scCulturePredict
#'
#' @description
#' Loads all required packages for the scCulturePredict package. If a package is not installed,
#' it will be installed automatically.
#'
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return Invisible NULL. This function is called for its side effects of loading the required packages into the current R session.
#' @export
#'
#' @examples
#' load_packages()
load_packages <- function(verbose = TRUE) {
  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  required <- c("Seurat", "dplyr", "ggplot2", "tidyverse", "MASS", "e1071", "caret")

  for (pkg in required) {
    tryCatch(
      {
        if (!pkg %in% rownames(installed.packages())) {
          warning(sprintf("Package %s is not installed. Please install it manually: install.packages('%s')", pkg, pkg))
        }
        if (verbose) message(sprintf("Loading package: %s", pkg))
        requireNamespace(pkg, quietly = TRUE)
      },
      error = function(e) {
        stop(sprintf("Failed to load package %s: %s", pkg, e$message))
      }
    )
  }
}

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

#' Prepare files for Seurat compatibility
#'
#' @description
#' Prepares input files for Seurat by copying and renaming them to match Seurat's requirements.
#' Also handles malformed files by removing problematic first rows.
#'
#' @examples
#' \donttest{
#' # Example with temporary directory
#' tmp_dir <- tempdir()
#' # Assume files exist in tmp_dir
#' # prepare_files_for_seurat(tmp_dir, "experiment1")
#' }
#' @param input_dir Character string specifying the directory containing the original files
#' @param output_dir Character string specifying the directory for prepared files
#' @param experiment_id Character string specifying the experiment ID prefix in filenames
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return Character string with the path to the prepared directory
#' @keywords internal
prepare_files_for_seurat <- function(input_dir, output_dir, experiment_id, verbose = TRUE) {
  # Input validation
  if (!is.character(input_dir) || length(input_dir) != 1) {
    stop("input_dir must be a single character string")
  }
  if (!dir.exists(input_dir)) {
    stop(sprintf("Input directory not found: %s", input_dir))
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }

  if (!is.character(experiment_id) || length(experiment_id) != 1) {
    stop("experiment_id must be a single character string")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  tryCatch(
    {
      # Create output directory if it doesn't exist
      if (!dir.exists(output_dir)) {
        if (verbose) message(sprintf("Creating output directory: %s", output_dir))
        dir.create(output_dir, recursive = TRUE)
      }

      # Define file mappings
      file_mappings <- list(
        barcodes = c(
          paste0(experiment_id, "_barcodes.tsv.gz"),
          "barcodes.tsv.gz"
        ),
        features = c(
          paste0(experiment_id, "_features.tsv.gz"),
          "features.tsv.gz"
        ),
        matrix = c(
          paste0(experiment_id, "_yeastdropseq_dge.mtx"),
          "matrix.mtx.gz"
        ),
        metadata = c(
          paste0(experiment_id, "_metadata.tsv.gz"),
          "metadata.tsv.gz"
        )
      )

      # Process each file
      for (file_type in names(file_mappings)) {
        input_file <- file.path(input_dir, file_mappings[[file_type]][1])
        output_file <- file.path(output_dir, file_mappings[[file_type]][2])

        if (file.exists(input_file)) {
          if (verbose) message(sprintf("Processing %s file...", file_type))

          # For barcodes and features, remove first row if it contains 'x'
          if (file_type %in% c("barcodes", "features")) {
            lines <- readLines(input_file, n = 5)
            if (grepl("^x\\s*$", lines[1]) || grepl("^x\\s+.*$", lines[1])) {
              if (verbose) message(sprintf("Detected malformed header in %s, fixing...", file_type))
              # Read file skipping first line
              data <- read.csv(input_file, sep = "\t", header = (file_type == "features"))
              # Write to new file
              write.table(data, gsub("\\.gz$", "", output_file),
                sep = "\t", row.names = FALSE, col.names = (file_type == "features")
              )
              # Compress (conditional on R.utils availability)
              if (requireNamespace("R.utils", quietly = TRUE)) {
                R.utils::gzip(gsub("\\.gz$", "", output_file), destname = output_file, remove = TRUE)
              } else {
                warning("R.utils package not available. File will not be compressed.")
                # Rename uncompressed file to match expected output
                file.rename(gsub("\\.gz$", "", output_file), output_file)
              }
            } else {
              # Just copy the file
              if (verbose) message(sprintf("Copying %s file...", file_type))
              file.copy(input_file, output_file)
            }
          } else {
            # For matrix and metadata, just copy
            if (verbose) message(sprintf("Copying %s file...", file_type))
            file.copy(input_file, output_file)
          }
        } else {
          warning(sprintf("File not found: %s", input_file))
        }
      }

      return(output_dir)
    },
    error = function(e) {
      stop(sprintf("Failed to prepare files for Seurat: %s", e$message))
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
#' @param use_shell_script Logical indicating whether to use shell script for file preparation.
#'   If TRUE, attempts to use shell script first, falling back to R-based preparation if it fails.
#'   Default is TRUE.
#' @param metadata_file Character string specifying the name of the metadata file.
#'   If NULL, will be automatically detected based on experiment_id.
#' @param counts_file Character string specifying the name of the counts file.
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
#' # Example with temporary files
#' # Create temporary directory and files
#' temp_dir <- tempdir()
#' exp_id <- "test_exp"
#'
#' # Create mock counts file
#' counts_data <- matrix(rpois(500, 5), nrow = 50)
#' rownames(counts_data) <- paste0("Gene", seq_len(50))
#' colnames(counts_data) <- paste0("Cell", seq_len(10))
#' counts_file <- file.path(temp_dir, paste0(exp_id, "_counts.csv"))
#' write.csv(counts_data, counts_file)
#'
#' # Create mock metadata file
#' metadata <- data.frame(
#'   cell_id = colnames(counts_data),
#'   sample = rep(c("Sample1", "Sample2"), each = 5)
#' )
#' metadata_file <- file.path(temp_dir, paste0(exp_id, "_metadata.csv"))
#' write.csv(metadata, metadata_file, row.names = FALSE)
#'
#' # Load the data
#' seurat_obj <- load_data(
#'   data_dir = temp_dir,
#'   experiment_id = exp_id,
#'   use_shell_script = FALSE
#' )
#'
#' # Clean up
#' unlink(counts_file)
#' unlink(metadata_file)
#' @seealso
#' \code{\link{preprocess_data}} for preprocessing the loaded data
#' \code{\link{reduce_dimensions}} for dimensionality reduction
#'
#' @export
load_data <- function(data_dir,
                      experiment_id,
                      use_shell_script = TRUE,
                      metadata_file = NULL,
                      counts_file = NULL,
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

  if (!is.null(counts_file) && (!is.character(counts_file) || length(counts_file) != 1)) {
    stop("counts_file must be a single character string or NULL")
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

  # Auto-detect data format (10X vs CSV)
  if (verbose) message("Detecting data format...")

  # Check for 10X format files
  tenx_files <- c(
    "matrix.mtx", "matrix.mtx.gz", "barcodes.tsv", "barcodes.tsv.gz",
    "features.tsv", "features.tsv.gz", "genes.tsv", "genes.tsv.gz"
  )
  tenx_present <- any(file.exists(file.path(data_dir, tenx_files)))

  if (tenx_present) {
    # 10X Genomics format detected
    if (verbose) message("10X Genomics format detected")

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

    data_format <- "10X"
  } else {
    # CSV format
    if (verbose) message("CSV format detected")

    # Find metadata file for CSV format
    if (is.null(metadata_file)) {
      metadata_file <- list.files(data_dir,
        pattern = paste0(experiment_id, ".*metadata.*\\.csv$"),
        full.names = TRUE
      )
      if (length(metadata_file) == 0) {
        stop(sprintf("No metadata file found for experiment %s", experiment_id))
      }
      if (length(metadata_file) > 1) {
        warning(sprintf("Multiple metadata files found, using: %s", metadata_file[1]))
        metadata_file <- metadata_file[1]
      }
    } else {
      metadata_file <- file.path(data_dir, metadata_file)
      if (!file.exists(metadata_file)) {
        stop(sprintf("Metadata file not found: %s", metadata_file))
      }
    }

    # Find counts file for CSV format
    if (is.null(counts_file)) {
      counts_file <- list.files(data_dir,
        pattern = paste0(experiment_id, ".*counts.*\\.csv$"),
        full.names = TRUE
      )
      if (length(counts_file) == 0) {
        stop(sprintf("No counts file found for experiment %s", experiment_id))
      }
      if (length(counts_file) > 1) {
        warning(sprintf("Multiple counts files found, using: %s", counts_file[1]))
        counts_file <- counts_file[1]
      }
    } else {
      counts_file <- file.path(data_dir, counts_file)
      if (!file.exists(counts_file)) {
        stop(sprintf("Counts file not found: %s", counts_file))
      }
    }

    data_format <- "CSV"
  }

  if (data_format == "10X") {
    # Load 10X Genomics format (following original working approach)
    if (verbose) message("Loading 10X Genomics data...")

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
  } else {
    # Load CSV format
    if (verbose) message("Loading CSV format data...")

    # Load and validate metadata
    metadata <- tryCatch(
      {
        read.csv(metadata_file, row.names = 1)
      },
      error = function(e) {
        stop(sprintf("Error loading metadata file: %s", e$message))
      }
    )

    if (!"sample" %in% colnames(metadata)) {
      stop("Metadata must contain a 'sample' column")
    }

    # Load and validate counts
    if (verbose) message("Loading counts data...")
    counts <- tryCatch(
      {
        read.csv(counts_file, row.names = 1)
      },
      error = function(e) {
        stop(sprintf("Error loading counts file: %s", e$message))
      }
    )

    # Validate data dimensions
    if (ncol(counts) != nrow(metadata)) {
      stop(sprintf(
        "Number of cells in counts (%d) does not match metadata (%d)",
        ncol(counts), nrow(metadata)
      ))
    }

    # Create Seurat object
    if (verbose) message("Creating Seurat object...")
    seurat_obj <- tryCatch(
      {
        Seurat::CreateSeuratObject(
          counts = counts,
          meta.data = metadata,
          min.cells = min_cells,
          min.features = min_features
        )
      },
      error = function(e) {
        stop(sprintf("Error creating Seurat object: %s", e$message))
      }
    )
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
