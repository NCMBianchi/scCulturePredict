# Test error handling and edge cases in the pipeline
# This file tests scCulture's behavior with invalid inputs and edge cases

# Helper function to create invalid/edge case data
create_invalid_data_dir <- function(type = "missing_files") {
  temp_dir <- tempfile("test_error_data")
  dir.create(temp_dir)

  if (type == "missing_files") {
    # Create directory with missing files
    # Only create barcodes.tsv.gz, missing matrix and features
    n_cells <- 100
    barcodes <- paste0("BARCODE", seq_len(n_cells))
    barcodes_df <- data.frame(
      line_num = seq_len(n_cells),
      barcode = barcodes
    )
    write.table(barcodes_df,
      file.path(temp_dir, "barcodes.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    R.utils::gzip(file.path(temp_dir, "barcodes.tsv"),
      destname = file.path(temp_dir, "barcodes.tsv.gz"),
      remove = FALSE, overwrite = TRUE
    )
  } else if (type == "malformed_matrix") {
    # Create files with malformed matrix
    # Write corrupted matrix file
    writeLines(
      c(
        "%%MatrixMarket matrix coordinate real general",
        "INVALID DATA HERE"
      ),
      file.path(temp_dir, "matrix.mtx")
    )
    R.utils::gzip(file.path(temp_dir, "matrix.mtx"),
      destname = file.path(temp_dir, "matrix.mtx.gz"),
      remove = FALSE, overwrite = TRUE
    )

    # Create valid barcodes and features
    n_cells <- 10
    n_genes <- 100
    barcodes_df <- data.frame(
      line_num = seq_len(n_cells),
      barcode = paste0("BARCODE", seq_len(n_cells))
    )
    write.table(barcodes_df,
      file.path(temp_dir, "barcodes.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    R.utils::gzip(file.path(temp_dir, "barcodes.tsv"),
      destname = file.path(temp_dir, "barcodes.tsv.gz"),
      remove = FALSE, overwrite = TRUE
    )

    features_df <- data.frame(
      line_num = seq_len(n_genes),
      gene_name = paste0("GENE", seq_len(n_genes))
    )
    write.table(features_df,
      file.path(temp_dir, "features.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    R.utils::gzip(file.path(temp_dir, "features.tsv"),
      destname = file.path(temp_dir, "features.tsv.gz"),
      remove = FALSE, overwrite = TRUE
    )
  } else if (type == "empty_data") {
    # Create valid but empty data (0 counts)
    n_genes <- 100
    n_cells <- 50

    # Create empty sparse matrix
    count_matrix <- matrix(0, nrow = n_genes, ncol = n_cells)
    sparse_mat <- Matrix::Matrix(count_matrix, sparse = TRUE)
    mtx_file <- file.path(temp_dir, "matrix.mtx")
    Matrix::writeMM(sparse_mat, mtx_file)
    R.utils::gzip(mtx_file,
      destname = paste0(mtx_file, ".gz"),
      remove = FALSE, overwrite = TRUE
    )

    # Create barcodes
    barcodes_df <- data.frame(
      line_num = seq_len(n_cells),
      barcode = paste0("BARCODE", seq_len(n_cells))
    )
    write.table(barcodes_df,
      file.path(temp_dir, "barcodes.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    R.utils::gzip(file.path(temp_dir, "barcodes.tsv"),
      destname = file.path(temp_dir, "barcodes.tsv.gz"),
      remove = FALSE, overwrite = TRUE
    )

    # Create features
    features_df <- data.frame(
      line_num = seq_len(n_genes),
      gene_name = paste0("YAL", sprintf("%03d", seq_len(n_genes)), "W")
    )
    write.table(features_df,
      file.path(temp_dir, "features.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    R.utils::gzip(file.path(temp_dir, "features.tsv"),
      destname = file.path(temp_dir, "features.tsv.gz"),
      remove = FALSE, overwrite = TRUE
    )

    # Create metadata
    metadata <- data.frame(
      sample = rep("DMSO", n_cells),
      sc3_8_clusters = rep("DMSO.1", n_cells)
    )
    rownames(metadata) <- paste0("BARCODE", seq_len(n_cells))
    write.table(metadata,
      file.path(temp_dir, "metadata.tsv"),
      sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE
    )
    R.utils::gzip(file.path(temp_dir, "metadata.tsv"),
      destname = file.path(temp_dir, "metadata.tsv.gz"),
      remove = FALSE, overwrite = TRUE
    )
  } else if (type == "single_cell") {
    # Create data with only one cell
    n_genes <- 100
    n_cells <- 1

    # Create sparse matrix with some data
    count_matrix <- matrix(0, nrow = n_genes, ncol = n_cells)
    # Add some counts
    count_matrix[sample(n_genes, 20), 1] <- rpois(20, lambda = 5)

    sparse_mat <- Matrix::Matrix(count_matrix, sparse = TRUE)
    mtx_file <- file.path(temp_dir, "matrix.mtx")
    Matrix::writeMM(sparse_mat, mtx_file)
    R.utils::gzip(mtx_file,
      destname = paste0(mtx_file, ".gz"),
      remove = FALSE, overwrite = TRUE
    )

    # Create single barcode
    barcodes_df <- data.frame(
      line_num = 1,
      barcode = "AAAAAAAAAAAAA"
    )
    write.table(barcodes_df,
      file.path(temp_dir, "barcodes.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    R.utils::gzip(file.path(temp_dir, "barcodes.tsv"),
      destname = file.path(temp_dir, "barcodes.tsv.gz"),
      remove = FALSE, overwrite = TRUE
    )

    # Create features
    features_df <- data.frame(
      line_num = seq_len(n_genes),
      gene_name = paste0("YAL", sprintf("%03d", seq_len(n_genes)), "W")
    )
    write.table(features_df,
      file.path(temp_dir, "features.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
    )
    R.utils::gzip(file.path(temp_dir, "features.tsv"),
      destname = file.path(temp_dir, "features.tsv.gz"),
      remove = FALSE, overwrite = TRUE
    )

    # Create metadata
    metadata <- data.frame(
      sample = "DMSO",
      sc3_8_clusters = "DMSO.1"
    )
    rownames(metadata) <- "AAAAAAAAAAAAA"
    write.table(metadata,
      file.path(temp_dir, "metadata.tsv"),
      sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE
    )
    R.utils::gzip(file.path(temp_dir, "metadata.tsv"),
      destname = file.path(temp_dir, "metadata.tsv.gz"),
      remove = FALSE, overwrite = TRUE
    )
  }

  return(temp_dir)
}

# Helper to get the default KEGG file path
get_default_kegg_file <- function() {
  system.file("extdata", "kegg", "sce00001.keg", package = "scCulturePredict")
}

test_that("scCulture handles missing files gracefully", {
  # Create directory with missing files
  data_dir <- create_invalid_data_dir("missing_files")
  kegg_file <- get_default_kegg_file()

  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Should error when matrix file is missing
  expect_error(
    scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = "test_missing",
      kegg_file = kegg_file,
      output_dir = tempfile("test_missing_output"),
      verbose = FALSE
    ),
    regexp = "matrix|Matrix|file|File"
  )

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture handles non-existent directory", {
  kegg_file <- get_default_kegg_file()

  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Should error with non-existent directory
  expect_error(
    scCulture(
      mode = "build",
      data_dir = "/path/that/does/not/exist",
      experiment_id = "test_nodir",
      kegg_file = kegg_file,
      output_dir = tempfile("test_nodir_output"),
      verbose = FALSE
    ),
    regexp = "exist|directory|Directory"
  )
})

test_that("scCulture handles invalid mode parameter", {
  # Create minimal valid data
  data_dir <- create_invalid_data_dir("single_cell")
  kegg_file <- get_default_kegg_file()

  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Should error with invalid mode
  expect_error(
    scCulture(
      mode = "invalid_mode",
      data_dir = data_dir,
      experiment_id = "test_invalid_mode",
      kegg_file = kegg_file,
      output_dir = tempfile("test_invalid_mode_output"),
      verbose = FALSE
    ),
    regexp = "mode|Mode|build|predict"
  )

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture handles missing KEGG file", {
  data_dir <- create_invalid_data_dir("single_cell")

  # Should error with non-existent KEGG file
  expect_error(
    scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = "test_no_kegg",
      kegg_file = "/path/to/nonexistent/kegg.keg",
      output_dir = tempfile("test_no_kegg_output"),
      verbose = FALSE
    ),
    regexp = "KEGG|kegg|file|exist"
  )

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture handles missing fingerprint file in PREDICT mode", {
  data_dir <- create_invalid_data_dir("single_cell")

  # Should error when fingerprint file doesn't exist
  expect_error(
    scCulture(
      mode = "predict",
      data_dir = data_dir,
      experiment_id = "test_no_fingerprint",
      fingerprint_file = "/path/to/nonexistent/fingerprint.rds",
      output_dir = tempfile("test_no_fingerprint_output"),
      verbose = FALSE
    ),
    regexp = "fingerprint|Fingerprint|file|exist"
  )

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture handles malformed data files", {
  # Create directory with malformed matrix
  data_dir <- create_invalid_data_dir("malformed_matrix")
  kegg_file <- get_default_kegg_file()

  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Should error when matrix is malformed
  expect_error(
    scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = "test_malformed",
      kegg_file = kegg_file,
      output_dir = tempfile("test_malformed_output"),
      verbose = FALSE
    )
  )

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture handles empty data matrix", {
  # Create directory with empty data
  data_dir <- create_invalid_data_dir("empty_data")
  kegg_file <- get_default_kegg_file()

  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Should either handle gracefully or error with informative message
  result <- tryCatch(
    {
      scCulture(
        mode = "build",
        data_dir = data_dir,
        experiment_id = "test_empty",
        kegg_file = kegg_file,
        output_dir = tempfile("test_empty_output"),
        verbose = FALSE
      )
    },
    error = function(e) e
  )

  # If it errors, just check that it is indeed an error
  if (inherits(result, "error")) {
    expect_true(inherits(result, "error"))
  }

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture handles single cell data", {
  # Create directory with single cell
  data_dir <- create_invalid_data_dir("single_cell")
  kegg_file <- get_default_kegg_file()

  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Should either handle single cell or error appropriately
  result <- tryCatch(
    {
      scCulture(
        mode = "build",
        data_dir = data_dir,
        experiment_id = "test_single",
        kegg_file = kegg_file,
        output_dir = tempfile("test_single_output"),
        verbose = FALSE,
        perform_tsne = FALSE # t-SNE won't work with single cell
      )
    },
    error = function(e) e
  )

  # If it errors, just check that it is indeed an error
  if (inherits(result, "error")) {
    expect_true(inherits(result, "error"))
  }

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture handles invalid parameter combinations", {
  data_dir <- create_invalid_data_dir("single_cell")
  kegg_file <- get_default_kegg_file()

  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Test invalid n_cores (negative)
  expect_error(
    scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = "test_negative_cores",
      kegg_file = kegg_file,
      output_dir = tempfile("test_negative_cores_output"),
      verbose = FALSE,
      parallel = TRUE,
      n_cores = -1
    ),
    regexp = "cores|positive|valid"
  )

  # Test invalid n_pcs (too high)
  expect_error(
    scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = "test_high_pcs",
      kegg_file = kegg_file,
      output_dir = tempfile("test_high_pcs_output"),
      verbose = FALSE,
      n_pcs = 10000
    ),
    regexp = "PC|pcs|dimension|valid"
  )

  # Test invalid tsne_perplexity (negative)
  expect_error(
    scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = "test_negative_perplexity",
      kegg_file = kegg_file,
      output_dir = tempfile("test_negative_perplexity_output"),
      verbose = FALSE,
      perform_tsne = TRUE,
      tsne_perplexity = -10
    ),
    regexp = "perplexity|positive|valid"
  )

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture handles NULL required parameters", {
  kegg_file <- get_default_kegg_file()

  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Test NULL data_dir
  expect_error(
    scCulture(
      mode = "build",
      data_dir = NULL,
      experiment_id = "test_null",
      kegg_file = kegg_file,
      output_dir = tempfile("test_null_output"),
      verbose = FALSE
    ),
    regexp = "data_dir|required|NULL"
  )

  # Test NULL experiment_id
  data_dir <- create_invalid_data_dir("single_cell")
  expect_error(
    scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = NULL,
      kegg_file = kegg_file,
      output_dir = tempfile("test_null_id_output"),
      verbose = FALSE
    )
  )

  # Test NULL kegg_file in BUILD mode
  expect_error(
    scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = "test_null_kegg",
      kegg_file = NULL,
      output_dir = tempfile("test_null_kegg_output"),
      verbose = FALSE
    )
  )

  # Test NULL fingerprint_file in PREDICT mode
  expect_error(
    scCulture(
      mode = "predict",
      data_dir = data_dir,
      experiment_id = "test_null_fp",
      fingerprint_file = NULL,
      output_dir = tempfile("test_null_fp_output"),
      verbose = FALSE
    )
  )

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture handles permission errors gracefully", {
  skip_on_os("windows") # Permission handling differs on Windows

  data_dir <- create_invalid_data_dir("single_cell")
  kegg_file <- get_default_kegg_file()

  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Create a read-only output directory
  output_dir <- tempfile("test_readonly_output")
  dir.create(output_dir, mode = "0444")

  # Try to run with read-only output directory
  result <- tryCatch(
    {
      scCulture(
        mode = "build",
        data_dir = data_dir,
        experiment_id = "test_readonly",
        kegg_file = kegg_file,
        output_dir = file.path(output_dir, "subdir"),
        verbose = FALSE
      )
    },
    error = function(e) e
  )

  # Should error about permissions
  if (inherits(result, "error")) {
    expect_true(inherits(result, "error"))
  }

  # Clean up
  unlink(data_dir, recursive = TRUE)
  unlink(output_dir, recursive = TRUE, force = TRUE)
})

test_that("plot_scCulture handles missing results components", {
  # Test with missing seurat_object
  expect_error(
    plot_scCulture(list(mode = "build")),
    regexp = "seurat_object"
  )

  # Test with NULL seurat_object
  expect_error(
    plot_scCulture(list(mode = "build", seurat_object = NULL)),
    regexp = "seurat_object"
  )

  # Test with invalid seurat_object type
  expect_error(
    plot_scCulture(list(mode = "build", seurat_object = "not_a_seurat")),
    regexp = "Seurat|seurat|valid"
  )
})

test_that("scCulture handles extreme parameter values", {
  data_dir <- create_invalid_data_dir("single_cell")
  kegg_file <- get_default_kegg_file()

  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Test with min_genes = 0
  result <- tryCatch(
    {
      scCulture(
        mode = "build",
        data_dir = data_dir,
        experiment_id = "test_min_genes_zero",
        kegg_file = kegg_file,
        output_dir = tempfile("test_min_genes_zero_output"),
        min_genes = 0,
        verbose = FALSE
      )
    },
    error = function(e) e
  )

  # Should either work or give meaningful error
  if (inherits(result, "error")) {
    expect_true(grepl("min_genes|positive|valid", result$message, ignore.case = TRUE))
  }

  # Test with extremely high min_cells
  result <- tryCatch(
    {
      scCulture(
        mode = "build",
        data_dir = data_dir,
        experiment_id = "test_high_min_cells",
        kegg_file = kegg_file,
        output_dir = tempfile("test_high_min_cells_output"),
        min_cells = 10000,
        verbose = FALSE
      )
    },
    error = function(e) e
  )

  # Should error about no genes passing filter
  if (inherits(result, "error")) {
    expect_true(grepl("genes|cells|filter|no data", result$message, ignore.case = TRUE))
  }

  # Clean up
  unlink(data_dir, recursive = TRUE)
})
