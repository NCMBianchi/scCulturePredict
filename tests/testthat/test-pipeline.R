# Test pipeline functions

# Create a mock data directory structure
create_mock_data <- function() {
  # Create temporary directory
  temp_dir <- tempfile("test_data")
  dir.create(temp_dir)
  # Directory cleanup handled by test framework

  # Create mock CSV format data (simpler and more reliable for testing)
  n_genes <- 200
  n_cells <- 100

  # Generate count matrix
  set.seed(123)
  count_matrix <- matrix(
    rpois(n_genes * n_cells, lambda = 5),
    nrow = n_genes, ncol = n_cells
  )

  # Create gene names that will match KEGG pathways
  gene_names <- c(
    paste0("YAL", sprintf("%03d", 1:50), "W"), # 50 genes
    paste0("YBR", sprintf("%03d", 1:50), "C"), # 50 genes
    paste0("YCL", sprintf("%03d", 1:50), "W"), # 50 genes
    paste0("YDR", sprintf("%03d", 1:50), "C") # 50 genes
  )

  # Create cell names
  cell_names <- paste0("CELL", 1:n_cells)

  # Add row and column names
  rownames(count_matrix) <- gene_names
  colnames(count_matrix) <- cell_names

  # Write counts file in CSV format
  counts_df <- as.data.frame(count_matrix)
  write.csv(counts_df, file.path(temp_dir, "test_experiment_counts.csv"))

  # Create metadata with proper structure
  metadata <- data.frame(
    row.names = cell_names,
    sample = rep(c("DMSO", "Guanine", "MPA", "Control"), length.out = n_cells),
    batch = rep(c("Batch1", "Batch2"), length.out = n_cells),
    nCount_RNA = colSums(count_matrix),
    nFeature_RNA = colSums(count_matrix > 0)
  )
  write.csv(metadata, file.path(temp_dir, "test_experiment_metadata.csv"))

  # Create mock KEGG file with matching genes
  writeLines(
    c(
      "C    Glycolysis / Gluconeogenesis",
      paste0("D    ", gene_names[1:10]),
      "C    Citrate cycle (TCA cycle)",
      paste0("D    ", gene_names[11:25]),
      "C    Pentose phosphate pathway",
      paste0("D    ", gene_names[26:40])
    ),
    file.path(temp_dir, "sce00001.keg")
  )

  return(temp_dir)
}

test_that("scCulture works correctly", {
  # Create mock data
  data_dir <- create_mock_data()

  # Test basic pipeline with explicit KEGG file
  results <- scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "test_experiment",
    kegg_file = file.path(data_dir, "sce00001.keg"),
    output_dir = tempfile("test_output"),
    use_shell_script = FALSE,
    perform_tsne = FALSE
  )

  # Check results structure
  expect_type(results, "list")
  expect_s4_class(results$seurat_object, "Seurat")
  expect_type(results$pathway_results, "list")
  expect_type(results$prediction_results, "list")
  expect_type(results$evaluation_results, "list")

  # Check Seurat object
  expect_true(all(c("UMAP_1", "UMAP_2") %in%
    names(results$seurat_object@meta.data)))
  expect_true(all(c(
    "predicted_sample_1", "predicted_sample_2",
    "classification_pred"
  ) %in%
    names(results$seurat_object@meta.data)))

  # Test with default KEGG file
  results_default <- scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "test_experiment",
    output_dir = tempfile("test_output_default"),
    use_shell_script = FALSE,
    perform_tsne = FALSE
  )

  # Check results structure with default KEGG
  expect_type(results_default, "list")
  expect_s4_class(results_default$seurat_object, "Seurat")
})

test_that("scCulture with progress works correctly", {
  # Create mock data
  data_dir <- create_mock_data()

  # Test basic pipeline with progress
  results <- scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "test_experiment",
    kegg_file = file.path(data_dir, "sce00001.keg"),
    output_dir = tempfile("test_output"),
    use_shell_script = FALSE,
    perform_tsne = FALSE,
    parallel = FALSE
  )

  # Check results structure
  expect_type(results, "list")
  expect_s4_class(results$seurat_object, "Seurat")
  expect_type(results$pathway_results, "list")
  expect_type(results$prediction_results, "list")
  expect_type(results$evaluation_results, "list")

  # Test parallel processing
  results_parallel <- scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "test_experiment",
    kegg_file = file.path(data_dir, "sce00001.keg"),
    output_dir = tempfile("test_output"),
    use_shell_script = FALSE,
    perform_tsne = FALSE,
    parallel = TRUE,
    n_cores = 2
  )

  # Check parallel results
  expect_type(results_parallel, "list")
  expect_s4_class(results_parallel$seurat_object, "Seurat")
})

test_that("scCulture error handling works correctly", {
  # Test with invalid data directory
  expect_error(scCulture(
    mode = "build",
    data_dir = "nonexistent_dir",
    experiment_id = "test_experiment",
    kegg_file = "sce00001.keg",
    output_dir = tempfile("test_output")
  ))

  # Test with invalid experiment ID (no longer needed as we accept any experiment_id)
  # Test with missing data files instead
  empty_dir <- tempfile("empty_data")
  dir.create(empty_dir)
  expect_error(scCulture(
    mode = "build",
    data_dir = empty_dir,
    experiment_id = "test_experiment",
    output_dir = tempfile("test_output")
  ))
  unlink(empty_dir, recursive = TRUE)

  # Test with invalid KEGG file
  expect_error(scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "test_experiment",
    kegg_file = "nonexistent.keg",
    output_dir = tempfile("test_output")
  ))
})
