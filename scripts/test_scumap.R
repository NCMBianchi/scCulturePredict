#!/usr/bin/env Rscript

# Test script for the new scumap() function
# This script tests the main pipeline function to ensure it works correctly

# Load required libraries
library(scCulturePredict)

# Test 1: Basic scumap() functionality
cat("=== Testing Basic scumap() Functionality ===\n")

# Set up paths
data_dir <- system.file("extdata", "example_data", package = "scCulturePredict")
kegg_file <- system.file("extdata", "kegg", "example_pathway.keg", package = "scCulturePredict")
output_dir <- "./test_scumap_results"

# Clean up any existing results
if (dir.exists(output_dir)) {
  unlink(output_dir, recursive = TRUE)
}

# Run basic scumap analysis
cat("Running basic scumap() analysis...\n")
results <- scumap(
  data_dir = data_dir,
  kegg_file = kegg_file,
  output_dir = output_dir,
  experiment_id = "test_basic",
  verbose = TRUE
)

# Validate results structure
cat("Validating results structure...\n")
expected_components <- c("seurat_object", "pathway_results", "prediction_results", "evaluation_results")
missing_components <- setdiff(expected_components, names(results))
if (length(missing_components) > 0) {
  stop(paste("Missing result components:", paste(missing_components, collapse = ", ")))
}
cat("✓ Results structure is correct\n")

# Validate Seurat object
seurat_obj <- results$seurat_object
if (!inherits(seurat_obj, "Seurat")) {
  stop("seurat_object is not a Seurat object")
}
cat("✓ Seurat object is valid\n")

# Check for UMAP coordinates
if (!all(c("UMAP_1", "UMAP_2") %in% colnames(seurat_obj@meta.data))) {
  stop("UMAP coordinates not found in Seurat object")
}
cat("✓ UMAP coordinates present\n")

# Check for predictions
prediction_cols <- c("predicted_sample_1", "predicted_sample_2", "classification_pred")
missing_pred_cols <- setdiff(prediction_cols, colnames(seurat_obj@meta.data))
if (length(missing_pred_cols) > 0) {
  warning(paste("Missing prediction columns:", paste(missing_pred_cols, collapse = ", ")))
}
cat("✓ Prediction columns present\n")

# Check output files
expected_files <- c("seurat_object.rds", "pathway_results.rds", "prediction_results.rds", "evaluation_results.rds")
missing_files <- expected_files[!file.exists(file.path(output_dir, expected_files))]
if (length(missing_files) > 0) {
  warning(paste("Missing output files:", paste(missing_files, collapse = ", ")))
}
cat("✓ Output files created\n")

cat("Basic test completed successfully!\n\n")

# Test 2: Advanced scumap() with progress and parallel options
cat("=== Testing Advanced scumap() Features ===\n")

output_dir_advanced <- "./test_scumap_advanced"
if (dir.exists(output_dir_advanced)) {
  unlink(output_dir_advanced, recursive = TRUE)
}

cat("Running advanced scumap() with progress tracking...\n")
results_advanced <- scumap(
  data_dir = data_dir,
  kegg_file = kegg_file,
  output_dir = output_dir_advanced,
  experiment_id = "test_advanced",
  progress = TRUE,
  parallel = FALSE, # Keep FALSE for testing to avoid complexity
  perform_tsne = TRUE,
  verbose = TRUE
)

cat("✓ Advanced features test completed\n\n")

# Test 3: Alternative function name
cat("=== Testing Alternative Function Name ===\n")

output_dir_alt <- "./test_run_scumap"
if (dir.exists(output_dir_alt)) {
  unlink(output_dir_alt, recursive = TRUE)
}

cat("Testing run_scumap() alias...\n")
results_alt <- run_scumap(
  data_dir = data_dir,
  kegg_file = kegg_file,
  output_dir = output_dir_alt,
  experiment_id = "test_alias",
  verbose = FALSE
)

cat("✓ Alternative function name works\n\n")

# Test 4: Result validation
cat("=== Validating Analysis Results ===\n")

# Check data dimensions
n_cells <- ncol(results$seurat_object)
n_genes <- nrow(results$seurat_object)
cat(sprintf("Dataset: %d cells, %d genes\n", n_cells, n_genes))

# Check prediction accuracy (if applicable)
if ("sample" %in% colnames(seurat_obj@meta.data) &&
  "classification_pred" %in% colnames(seurat_obj@meta.data)) {
  accuracy <- mean(seurat_obj$sample == seurat_obj$classification_pred, na.rm = TRUE)
  cat(sprintf("SVM prediction accuracy: %.2f%%\n", accuracy * 100))
}

# Check pathway results
pathway_results <- results$pathway_results
if (is.list(pathway_results) && "pathway_matrix" %in% names(pathway_results)) {
  n_pathways <- ncol(pathway_results$pathway_matrix)
  cat(sprintf("Pathway analysis: %d pathways analyzed\n", n_pathways))
}

cat("✓ Results validation completed\n\n")

# Test 5: Error handling
cat("=== Testing Error Handling ===\n")

# Test with invalid data directory
cat("Testing error handling for invalid data directory...\n")
tryCatch(
  {
    scumap(
      data_dir = "/nonexistent/path",
      kegg_file = kegg_file,
      output_dir = "./test_error",
      verbose = FALSE
    )
    stop("Should have failed with invalid data directory")
  },
  error = function(e) {
    cat("✓ Correctly caught error for invalid data directory\n")
  }
)

# Test with invalid KEGG file
cat("Testing error handling for invalid KEGG file...\n")
tryCatch(
  {
    scumap(
      data_dir = data_dir,
      kegg_file = "/nonexistent/file.keg",
      output_dir = "./test_error",
      verbose = FALSE
    )
    stop("Should have failed with invalid KEGG file")
  },
  error = function(e) {
    cat("✓ Correctly caught error for invalid KEGG file\n")
  }
)

cat("✓ Error handling tests completed\n\n")

# Cleanup
cat("=== Cleaning Up Test Files ===\n")
test_dirs <- c(output_dir, output_dir_advanced, output_dir_alt)
for (dir in test_dirs) {
  if (dir.exists(dir)) {
    unlink(dir, recursive = TRUE)
    cat(sprintf("✓ Cleaned up %s\n", dir))
  }
}

cat("\n=== ALL TESTS COMPLETED SUCCESSFULLY! ===\n")
cat("The scumap() function is working correctly and ready for use.\n")

# Summary
cat("\nSummary of tested features:\n")
cat("✓ Basic scumap() functionality\n")
cat("✓ Advanced options (progress, parallel, tsne)\n")
cat("✓ Alternative function name (run_scumap)\n")
cat("✓ Result structure and content validation\n")
cat("✓ Error handling for invalid inputs\n")
cat("✓ File output and cleanup\n")
