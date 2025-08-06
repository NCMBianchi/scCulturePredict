# Test scCulture pipeline function

# Helper function to create mock data matching real 10X format
create_mock_10x_data <- function() {
  # Create temporary directory
  temp_dir <- tempfile("test_data")
  dir.create(temp_dir)

  # Set dimensions that won't get filtered out
  n_genes <- 1000
  n_cells <- 500 # Enough cells to survive filtering

  # Generate sparse count matrix (~80% zeros like real scRNA-seq)
  set.seed(123)
  count_matrix <- matrix(0, nrow = n_genes, ncol = n_cells)

  # Fill ~20% of entries with count data
  n_nonzero <- round(n_genes * n_cells * 0.2)
  for (i in seq_len(n_nonzero)) {
    row_idx <- sample(n_genes, 1)
    col_idx <- sample(n_cells, 1)
    count_matrix[row_idx, col_idx] <- count_matrix[row_idx, col_idx] + rpois(1, lambda = 3)
  }

  # Ensure each cell has sufficient expressed genes (100-500)
  for (i in seq_len(n_cells)) {
    n_expressed <- sample(100:500, 1)
    expressed_genes <- sample(n_genes, n_expressed)
    for (gene in expressed_genes) {
      if (count_matrix[gene, i] == 0) {
        count_matrix[gene, i] <- rpois(1, lambda = 2) + 1
      }
    }
  }

  # Create realistic yeast gene names
  gene_names <- character(n_genes)
  gene_types <- c(
    "YAL", "YBR", "YCL", "YDR", "YEL", "YFR", "YGL", "YGR",
    "YHL", "YHR", "YIL", "YJL", "YKL", "YKR", "YLL", "YLR",
    "YML", "YMR", "YNL", "YNR", "YOL", "YOR", "YPL", "YPR"
  )

  for (i in seq_len(n_genes)) {
    prefix <- sample(gene_types, 1)
    number <- sprintf("%03d", sample(1:999, 1))
    suffix <- sample(c("W", "C"), 1)
    gene_names[i] <- paste0(prefix, number, suffix)
  }

  # Create realistic barcodes
  bases <- c("A", "C", "G", "T")
  barcodes <- character(n_cells)
  for (i in seq_len(n_cells)) {
    barcodes[i] <- paste0(sample(bases, 12, replace = TRUE), collapse = "")
  }

  # Write matrix.mtx in Matrix Market format
  sparse_mat <- Matrix::Matrix(count_matrix, sparse = TRUE)
  mtx_file <- file.path(temp_dir, "matrix.mtx")
  Matrix::writeMM(sparse_mat, mtx_file)

  # Also create gzipped version
  R.utils::gzip(mtx_file, destname = paste0(mtx_file, ".gz"), remove = FALSE, overwrite = TRUE)

  # Write barcodes.tsv with line numbers (matching real format)
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

  # Write features.tsv with line numbers (matching real format)
  features_df <- data.frame(
    line_num = seq_len(n_genes),
    gene_name = gene_names
  )
  write.table(features_df,
    file.path(temp_dir, "features.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE
  )
  R.utils::gzip(file.path(temp_dir, "features.tsv"),
    destname = file.path(temp_dir, "features.tsv.gz"),
    remove = FALSE, overwrite = TRUE
  )

  # Write metadata.tsv with header (matching real format)
  metadata <- data.frame(
    sample = rep(c("DMSO", "Guanine", "Adenine", "Thymine"), length.out = n_cells),
    sc3_8_clusters = rep(c("DMSO.1", "Guanine.1", "Adenine.1", "Thymine.1"), length.out = n_cells),
    sc3_4_clusters = rep(c("1", "2", "3", "4"), length.out = n_cells),
    PC_1 = rnorm(n_cells),
    PC_2 = rnorm(n_cells),
    PC_3 = rnorm(n_cells),
    UMAP_1 = rnorm(n_cells, sd = 5),
    UMAP_2 = rnorm(n_cells, sd = 5)
  )
  rownames(metadata) <- barcodes
  write.table(metadata,
    file.path(temp_dir, "metadata.tsv"),
    sep = "\t", quote = FALSE, row.names = TRUE, col.names = TRUE
  )
  R.utils::gzip(file.path(temp_dir, "metadata.tsv"),
    destname = file.path(temp_dir, "metadata.tsv.gz"),
    remove = FALSE, overwrite = TRUE
  )

  return(temp_dir)
}

# Helper to get the default KEGG file path
get_default_kegg_file <- function() {
  system.file("extdata", "kegg", "sce00001.keg", package = "scCulturePredict")
}

test_that("scCulture BUILD mode works with 10X format data", {
  # Create mock 10X data
  data_dir <- create_mock_10x_data()

  # Get the actual KEGG file from the package
  kegg_file <- get_default_kegg_file()

  # Skip test if KEGG file not found
  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Run scCulture in BUILD mode
  expect_no_error({
    results <- scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = "test_experiment",
      kegg_file = kegg_file, # Use real KEGG file
      output_dir = tempfile("test_output"),
      use_shell_script = FALSE,
      perform_tsne = FALSE,
      verbose = FALSE,
      progress = FALSE
    )
  })

  # Basic checks on the results
  expect_type(results, "list")
  expect_s4_class(results$seurat_object, "Seurat")
  expect_true("pathway_results" %in% names(results))
  expect_true("prediction_results" %in% names(results))
  expect_true("evaluation_results" %in% names(results))
  expect_true("fingerprint_file" %in% names(results))

  # Check that fingerprint file was created
  expect_true(file.exists(results$fingerprint_file))

  # Check Seurat object has necessary metadata
  seurat_meta <- results$seurat_object@meta.data
  expect_true("sample" %in% colnames(seurat_meta))
  expect_true("UMAP_1" %in% colnames(seurat_meta))
  expect_true("UMAP_2" %in% colnames(seurat_meta))

  # Check that enough cells survived filtering
  expect_gt(ncol(results$seurat_object), 50) # At least 50 cells should survive

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture PREDICT mode works with fingerprint file", {
  # Create mock 10X data for both BUILD and PREDICT
  build_dir <- create_mock_10x_data()
  predict_dir <- create_mock_10x_data()

  # Get the actual KEGG file
  kegg_file <- get_default_kegg_file()

  # Skip test if KEGG file not found
  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # First run BUILD mode to generate fingerprints
  build_results <- scCulture(
    mode = "build",
    data_dir = build_dir,
    experiment_id = "build_experiment",
    kegg_file = kegg_file,
    output_dir = tempfile("build_output"),
    use_shell_script = FALSE,
    perform_tsne = FALSE,
    verbose = FALSE,
    progress = FALSE
  )

  # Now run PREDICT mode using the fingerprints
  expect_no_error({
    predict_results <- scCulture(
      mode = "predict",
      data_dir = predict_dir,
      experiment_id = "predict_experiment",
      fingerprint_file = build_results$fingerprint_file,
      output_dir = tempfile("predict_output"),
      use_shell_script = FALSE,
      perform_tsne = FALSE,
      verbose = FALSE,
      progress = FALSE
    )
  })

  # Check PREDICT mode results
  expect_type(predict_results, "list")
  expect_s4_class(predict_results$seurat_object, "Seurat")
  expect_true("prediction_results" %in% names(predict_results))

  # Check predictions were added to metadata
  predict_meta <- predict_results$seurat_object@meta.data
  expect_true("predicted_sample_1" %in% colnames(predict_meta) ||
    "classification_pred" %in% colnames(predict_meta))

  # Clean up
  unlink(build_dir, recursive = TRUE)
  unlink(predict_dir, recursive = TRUE)
})

test_that("scCulture handles invalid inputs gracefully", {
  # Test with non-existent directory
  expect_error(
    scCulture(
      mode = "build",
      data_dir = "nonexistent_dir",
      experiment_id = "test",
      output_dir = tempfile("test_output")
    ),
    regexp = "not found|does not exist"
  )

  # Test with invalid mode
  data_dir <- create_mock_10x_data()
  expect_error(
    scCulture(
      mode = "invalid_mode",
      data_dir = data_dir,
      experiment_id = "test",
      output_dir = tempfile("test_output")
    ),
    regexp = "mode.*build.*predict"
  )

  # Test PREDICT mode without fingerprint file
  expect_error(
    scCulture(
      mode = "predict",
      data_dir = data_dir,
      experiment_id = "test",
      output_dir = tempfile("test_output")
    ),
    regexp = "fingerprint"
  )

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture works with default KEGG file", {
  # Create mock 10X data
  data_dir <- create_mock_10x_data()

  # Run without specifying KEGG file (should use default)
  expect_no_error({
    results <- scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = "test_experiment",
      output_dir = tempfile("test_output"),
      use_shell_script = FALSE,
      perform_tsne = FALSE,
      verbose = FALSE,
      progress = FALSE
    )
  })

  # Check that analysis completed
  expect_type(results, "list")
  expect_s4_class(results$seurat_object, "Seurat")

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("scCulture parallel processing works", {
  # Create mock 10X data
  data_dir <- create_mock_10x_data()

  # Get the actual KEGG file
  kegg_file <- get_default_kegg_file()

  # Skip test if KEGG file not found or if not enough cores
  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  if (parallel::detectCores() < 2) {
    skip("Not enough cores for parallel processing test")
  }

  # Run with parallel processing
  expect_no_error({
    results <- scCulture(
      mode = "build",
      data_dir = data_dir,
      experiment_id = "test_parallel",
      kegg_file = kegg_file,
      output_dir = tempfile("test_parallel_output"),
      use_shell_script = FALSE,
      perform_tsne = FALSE,
      verbose = FALSE,
      progress = FALSE,
      parallel = TRUE,
      n_cores = 2
    )
  })

  # Check results
  expect_type(results, "list")
  expect_s4_class(results$seurat_object, "Seurat")

  # Clean up
  unlink(data_dir, recursive = TRUE)
})
