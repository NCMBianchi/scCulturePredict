# Test pipeline functions

# Create a mock data directory structure
create_mock_data <- function() {
  # Create temporary directory
  temp_dir <- tempfile("test_data")
  dir.create(temp_dir)
  # Directory cleanup handled by test framework

  # Create mock data files
  count_matrix <- matrix(rpois(100, 10), nrow = 10, ncol = 10)
  rownames(count_matrix) <- paste0("gene", 1:10)
  colnames(count_matrix) <- paste0("cell", 1:10)
  write.csv(
    count_matrix,
    file.path(temp_dir, "GSE165686_counts.csv"),
    row.names = TRUE
  )

  # Create metadata with matching row names
  metadata <- data.frame(
    sample = rep(c("A", "B"), each = 5),
    row.names = paste0("cell", 1:10)
  )
  write.csv(
    metadata,
    file.path(temp_dir, "GSE165686_metadata.csv"),
    row.names = TRUE
  )

  # Create mock KEGG file
  writeLines(
    c(
      "A\tPATHWAY1",
      "B\tgene1",
      "B\tgene2",
      "A\tPATHWAY2",
      "B\tgene3",
      "B\tgene4"
    ),
    file.path(temp_dir, "sce00001.keg")
  )

  return(temp_dir)
}

test_that("scCulture works correctly", {
  # Create mock data
  data_dir <- create_mock_data()

  # Test basic pipeline
  results <- scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "GSE165686",
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
})

test_that("scCulture with progress works correctly", {
  # Create mock data
  data_dir <- create_mock_data()

  # Test basic pipeline with progress
  results <- scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "GSE165686",
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
    experiment_id = "GSE165686",
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
    experiment_id = "GSE165686",
    kegg_file = "sce00001.keg",
    output_dir = tempfile("test_output")
  ))

  # Test with invalid experiment ID
  data_dir <- create_mock_data()
  expect_error(scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "invalid_id",
    kegg_file = file.path(data_dir, "sce00001.keg"),
    output_dir = tempfile("test_output")
  ))

  # Test with invalid KEGG file
  expect_error(scCulture(
    mode = "build",
    data_dir = data_dir,
    experiment_id = "GSE165686",
    kegg_file = "nonexistent.keg",
    output_dir = tempfile("test_output")
  ))
})
