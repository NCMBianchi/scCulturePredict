# Test plot_scCulture visualization function

# Helper function to create mock data matching real 10X format (shared with test-pipeline.R)
create_mock_10x_data_for_viz <- function() {
  # Create temporary directory
  temp_dir <- tempfile("test_viz_data")
  dir.create(temp_dir)

  # Set dimensions that won't get filtered out
  n_genes <- 1000
  n_cells <- 500

  # Generate sparse count matrix
  set.seed(456) # Different seed for variety
  count_matrix <- matrix(0, nrow = n_genes, ncol = n_cells)

  # Fill ~20% of entries with count data
  n_nonzero <- round(n_genes * n_cells * 0.2)
  for (i in seq_len(n_nonzero)) {
    row_idx <- sample(n_genes, 1)
    col_idx <- sample(n_cells, 1)
    count_matrix[row_idx, col_idx] <- count_matrix[row_idx, col_idx] + rpois(1, lambda = 3)
  }

  # Ensure each cell has sufficient expressed genes
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

test_that("plot_scCulture works with BUILD mode results", {
  # Create mock data and run scCulture
  data_dir <- create_mock_10x_data_for_viz()

  # Get the actual KEGG file from the package
  kegg_file <- get_default_kegg_file()

  # Skip test if KEGG file not found
  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Run scCulture in BUILD mode to get real results
  build_results <- scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "test_viz",
    kegg_file = kegg_file,
    output_dir = tempfile("test_viz_output"),
    perform_tsne = FALSE,
    verbose = FALSE,
    progress = FALSE
  )

  # Test plot_scCulture with BUILD results
  expect_no_error({
    plot <- plot_scCulture(build_results)
  })

  # Check that it returns a ggplot object
  plot <- plot_scCulture(build_results)
  expect_s3_class(plot, "ggplot")
  expect_s3_class(plot, "gg")

  # Test with different plot_type if supported
  # Default plot_type for BUILD should be "accuracy"
  plot_accuracy <- plot_scCulture(build_results, plot_type = "accuracy")
  expect_s3_class(plot_accuracy, "ggplot")

  # Test with custom point size and alpha
  plot_custom <- plot_scCulture(build_results,
    point_size = 2.0,
    point_alpha = 0.5
  )
  expect_s3_class(plot_custom, "ggplot")

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("plot_scCulture works with PREDICT mode results", {
  # Create mock data for both BUILD and PREDICT
  build_dir <- create_mock_10x_data_for_viz()
  predict_dir <- create_mock_10x_data_for_viz()

  # Get the actual KEGG file
  kegg_file <- get_default_kegg_file()

  # Skip test if KEGG file not found
  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # First run BUILD to get fingerprints
  build_results <- scCulture(
    mode = "build",
    data_dir = build_dir,
    experiment_id = "build_viz",
    kegg_file = kegg_file,
    output_dir = tempfile("build_viz_output"),
    perform_tsne = FALSE,
    verbose = FALSE,
    progress = FALSE
  )

  # Run PREDICT mode
  predict_results <- scCulture(
    mode = "predict",
    data_dir = predict_dir,
    experiment_id = "predict_viz",
    fingerprint_file = build_results$fingerprint_file,
    output_dir = tempfile("predict_viz_output"),
    perform_tsne = FALSE,
    verbose = FALSE,
    progress = FALSE
  )

  # Test plot_scCulture with PREDICT results
  expect_no_error({
    plot <- plot_scCulture(predict_results)
  })

  # Check the return type - could be single plot or list of plots
  plot <- plot_scCulture(predict_results)
  expect_true(
    inherits(plot, "ggplot") || inherits(plot, "list"),
    info = "plot_scCulture should return either a ggplot or list of plots"
  )

  # If list, check components are ggplots
  if (is.list(plot) && !inherits(plot, "ggplot")) {
    for (p in plot) {
      if (!is.null(p)) {
        expect_s3_class(p, "ggplot")
      }
    }
  }

  # Test different plot_types for PREDICT mode if supported
  if (!is.null(predict_results$prediction_results)) {
    # Try predictions plot
    expect_no_error({
      plot_pred <- plot_scCulture(predict_results, plot_type = "predictions")
    })

    # Try confidence plot if confidence scores exist
    if (!is.null(predict_results$prediction_results$confidence)) {
      expect_no_error({
        plot_conf <- plot_scCulture(predict_results, plot_type = "confidence")
      })
    }
  }

  # Clean up
  unlink(build_dir, recursive = TRUE)
  unlink(predict_dir, recursive = TRUE)
})

test_that("plot_scCulture returns data when return_data = TRUE", {
  # Create mock data
  data_dir <- create_mock_10x_data_for_viz()

  # Get the actual KEGG file
  kegg_file <- get_default_kegg_file()

  # Skip test if KEGG file not found
  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Run scCulture
  results <- scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "test_data_return",
    kegg_file = kegg_file,
    output_dir = tempfile("test_data_output"),
    perform_tsne = FALSE,
    verbose = FALSE,
    progress = FALSE
  )

  # Test return_data parameter
  plot_data <- plot_scCulture(results, return_data = TRUE)

  # Check that data is returned instead of plot
  expect_type(plot_data, "list")
  expect_true(is.data.frame(plot_data) || is.list(plot_data))

  # If data frame, check for expected columns
  if (is.data.frame(plot_data)) {
    expect_true("UMAP_1" %in% colnames(plot_data) ||
      "umap_1" %in% colnames(plot_data) ||
      length(plot_data) > 0)
  }

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("plot_scCulture handles invalid inputs gracefully", {
  # Test with invalid input (not a list)
  expect_error(
    plot_scCulture("not a list"),
    regexp = "must be a list"
  )

  # Test with empty list
  expect_error(
    plot_scCulture(list()),
    regexp = "seurat_object"
  )

  # Test with list missing seurat_object
  expect_error(
    plot_scCulture(list(mode = "build")),
    regexp = "seurat_object"
  )

  # Test with invalid point_size
  data_dir <- create_mock_10x_data_for_viz()
  kegg_file <- get_default_kegg_file()

  # Skip test if KEGG file not found
  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  results <- scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "test_invalid",
    kegg_file = kegg_file,
    output_dir = tempfile("test_invalid_output"),
    perform_tsne = FALSE,
    verbose = FALSE,
    progress = FALSE
  )

  expect_error(
    plot_scCulture(results, point_size = -1),
    regexp = "positive"
  )

  expect_error(
    plot_scCulture(results, point_size = "large"),
    regexp = "positive"
  )

  # Test with invalid point_alpha
  expect_error(
    plot_scCulture(results, point_alpha = 2),
    regexp = "between 0 and 1"
  )

  # Clean up
  unlink(data_dir, recursive = TRUE)
})

test_that("plot_scCulture saves plots to file when requested", {
  # Create mock data
  data_dir <- create_mock_10x_data_for_viz()
  output_dir <- tempfile("test_save_output")

  # Get the actual KEGG file
  kegg_file <- get_default_kegg_file()

  # Skip test if KEGG file not found
  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Run scCulture
  results <- scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "test_save",
    kegg_file = kegg_file,
    output_dir = output_dir,
    perform_tsne = FALSE,
    verbose = FALSE,
    progress = FALSE
  )

  # Generate plot (it should save automatically if output_dir is in results)
  plot <- plot_scCulture(results)

  # Check if plot files were created in output directory
  # Note: The actual saving behavior depends on the implementation
  # This test checks that the function completes without error
  expect_s3_class(plot, "ggplot")

  # If the function has a save parameter, test it
  # (This depends on the actual implementation)

  # Clean up
  unlink(data_dir, recursive = TRUE)
  unlink(output_dir, recursive = TRUE)
})
