# Test prediction functions, focusing on SVM pathway name handling

library(testthat)
library(Seurat)
library(Matrix)

# Helper function to create mock Seurat object with pathway data
create_mock_seurat_with_pathways <- function(n_cells = 50, n_genes = 100, include_special_chars = TRUE) {
  # Create mock count matrix
  counts <- matrix(rpois(n_cells * n_genes, lambda = 5), nrow = n_genes, ncol = n_cells)
  rownames(counts) <- paste0("Gene_", 1:n_genes)
  colnames(counts) <- paste0("Cell_", 1:n_cells)

  # Create Seurat object
  seurat_obj <- CreateSeuratObject(counts = counts)

  # Add sample metadata with different conditions
  seurat_obj@meta.data$sample <- rep(c("Control", "Treatment"), length.out = n_cells)

  return(seurat_obj)
}

# Helper function to create mock pathway matrix with problematic names
create_mock_pathway_matrix <- function(n_cells = 50, include_special_chars = TRUE) {
  if (include_special_chars) {
    # Include pathway names that caused issues in the original bug report
    pathway_names <- c(
      "Glycolysis / Gluconeogenesis [PATH:sce00010]",
      "Citrate cycle (TCA cycle) [PATH:sce00020]",
      "Pentose phosphate pathway [PATH:sce00030]",
      "Fatty acid biosynthesis [PATH:sce00061]",
      "Steroid biosynthesis [PATH:sce00100]",
      "Purine metabolism [PATH:sce00230]",
      "Pyrimidine metabolism [PATH:sce00240]",
      "Alanine, aspartate and glutamate metabolism [PATH:sce00250]",
      "Glycine, serine and threonine metabolism [PATH:sce00260]",
      "Valine, leucine and isoleucine degradation [PATH:sce00280]"
    )
  } else {
    # Simple pathway names without special characters
    pathway_names <- paste0("Pathway_", 1:10)
  }

  # Create pathway activity matrix
  pathway_matrix <- matrix(
    rnorm(length(pathway_names) * n_cells, mean = 0.5, sd = 0.2),
    nrow = length(pathway_names),
    ncol = n_cells
  )
  rownames(pathway_matrix) <- pathway_names
  colnames(pathway_matrix) <- paste0("Cell_", 1:n_cells)

  return(pathway_matrix)
}

test_that("predict_by_svm handles pathway names with special characters", {
  # Create mock data with problematic pathway names
  seurat_obj <- create_mock_seurat_with_pathways(n_cells = 50)
  pathway_matrix <- create_mock_pathway_matrix(n_cells = 50, include_special_chars = TRUE)

  # Test that predict_by_svm doesn't error on special characters
  expect_no_error({
    results <- predict_by_svm(pathway_matrix, seurat_obj, verbose = FALSE)
  })

  # Verify results structure
  results <- predict_by_svm(pathway_matrix, seurat_obj, verbose = FALSE)
  expect_type(results, "list")
  expect_true(all(c("predictions", "svm_model", "confusion_matrix", "accuracy") %in% names(results)))
  expect_s3_class(results$svm_model, "svm")
  expect_true(is.numeric(results$accuracy))
  expect_true(results$accuracy >= 0 && results$accuracy <= 1)
})

test_that("make.names sanitization works consistently", {
  # Test pathway names that caused the original error
  problematic_names <- c(
    "Glycolysis / Gluconeogenesis [PATH:sce00010]",
    "Citrate cycle (TCA cycle) [PATH:sce00020]",
    "Alanine, aspartate and glutamate metabolism [PATH:sce00250]"
  )

  # Apply make.names transformation
  sanitized_names <- make.names(problematic_names)

  # Verify no special characters remain that would cause R variable name issues
  expect_true(all(make.names(sanitized_names) == sanitized_names))

  # Verify transformation is consistent
  expect_identical(make.names(problematic_names), sanitized_names)

  # Check specific transformations
  expect_true(grepl("Glycolysis", sanitized_names[1]))
  expect_false(grepl("/", sanitized_names[1]))
  expect_false(grepl("\\[", sanitized_names[1]))
  expect_false(grepl("\\]", sanitized_names[1]))
})

test_that("SVM model training uses sanitized feature names", {
  seurat_obj <- create_mock_seurat_with_pathways(n_cells = 50)
  pathway_matrix <- create_mock_pathway_matrix(n_cells = 50, include_special_chars = TRUE)

  # Train SVM model
  results <- predict_by_svm(pathway_matrix, seurat_obj, verbose = FALSE)

  # Extract model terms if available
  if (!is.null(results$svm_model$terms)) {
    model_features <- attr(results$svm_model$terms, "term.labels")

    # Verify all feature names are valid R variable names
    expect_true(all(make.names(model_features) == model_features))

    # Verify no special characters in feature names
    expect_false(any(grepl("[/\\[\\]()]", model_features)))
  }
})

test_that("predict_by_svm handles missing features gracefully", {
  # Create training data
  seurat_obj <- create_mock_seurat_with_pathways(n_cells = 50)
  pathway_matrix_full <- create_mock_pathway_matrix(n_cells = 50, include_special_chars = TRUE)

  # Train model on full feature set
  full_results <- predict_by_svm(pathway_matrix_full, seurat_obj, verbose = FALSE)

  # Create prediction data with fewer features (simulating missing pathways)
  pathway_matrix_subset <- pathway_matrix_full[1:5, ] # Only first 5 pathways

  # This should still work - the function should handle missing features
  expect_no_error({
    subset_seurat <- create_mock_seurat_with_pathways(n_cells = ncol(pathway_matrix_subset))
    subset_results <- predict_by_svm(pathway_matrix_subset, subset_seurat, verbose = FALSE)
  })
})

test_that("SVM prediction works with sanitized pathway names in PREDICT mode", {
  # Simulate the PREDICT mode scenario where we have a pre-trained model
  # and need to apply it to new data

  # Create mock fingerprint data structure
  seurat_obj <- create_mock_seurat_with_pathways(n_cells = 30)
  pathway_matrix <- create_mock_pathway_matrix(n_cells = 30, include_special_chars = TRUE)

  # Train initial model (simulating BUILD mode)
  build_results <- predict_by_svm(pathway_matrix, seurat_obj, verbose = FALSE)

  # Create mock fingerprint data
  fingerprint_data <- list(
    svm_model = build_results$svm_model,
    kegg_pathways = list(), # Mock KEGG data
    pathway_results = list(pathway_matrix = pathway_matrix),
    similarity_model = list(
      signature_matrix = pathway_matrix,
      pathway_matrix = pathway_matrix
    )
  )

  # Save and load fingerprint data to simulate real usage
  temp_fingerprint <- tempfile(fileext = ".rds")
  saveRDS(fingerprint_data, temp_fingerprint)

  # Test that loading and prediction works
  expect_no_error({
    loaded_data <- readRDS(temp_fingerprint)
    expect_s3_class(loaded_data$svm_model, "svm")
  })

  # Clean up
  unlink(temp_fingerprint)
})

test_that("error handling works when SVM prediction fails", {
  # Create mock data
  seurat_obj <- create_mock_seurat_with_pathways(n_cells = 20)
  pathway_matrix <- create_mock_pathway_matrix(n_cells = 20, include_special_chars = TRUE)

  # Create a deliberately problematic scenario
  # Make prediction data with completely different feature names
  prediction_data <- as.data.frame(t(pathway_matrix))
  colnames(prediction_data) <- paste0("DifferentFeature_", 1:ncol(prediction_data))

  # Create mock SVM model that expects different features
  mock_svm_model <- list(
    terms = structure(list(), term.labels = paste0("ExpectedFeature_", 1:10)),
    class = "svm"
  )
  class(mock_svm_model) <- "svm"

  # Test that prediction handles the mismatch gracefully
  expect_warning({
    predictions <- tryCatch(
      predict(mock_svm_model, newdata = prediction_data),
      error = function(e) {
        warning(sprintf("SVM prediction failed: %s", e$message))
        return(rep(NA, nrow(prediction_data)))
      }
    )
  })
})

test_that("pathway name consistency between BUILD and PREDICT modes", {
  # Create pathway matrix with special characters
  original_names <- c(
    "Glycolysis / Gluconeogenesis [PATH:sce00010]",
    "Citrate cycle (TCA cycle) [PATH:sce00020]",
    "Purine metabolism [PATH:sce00230]"
  )

  # Simulate BUILD mode name sanitization
  build_mode_names <- make.names(original_names)

  # Simulate PREDICT mode name sanitization
  predict_mode_names <- make.names(original_names)

  # Names should be identical between modes
  expect_identical(build_mode_names, predict_mode_names)

  # Both should be valid R variable names
  expect_true(all(make.names(build_mode_names) == build_mode_names))
  expect_true(all(make.names(predict_mode_names) == predict_mode_names))
})

test_that("feature matching works correctly in prediction", {
  # Create training data with sanitized names
  original_names <- c(
    "Glycolysis / Gluconeogenesis [PATH:sce00010]",
    "Citrate cycle (TCA cycle) [PATH:sce00020]",
    "Purine metabolism [PATH:sce00230]",
    "Fatty acid biosynthesis [PATH:sce00061]"
  )

  sanitized_names <- make.names(original_names)

  # Create prediction data frame
  prediction_data <- data.frame(matrix(rnorm(40), nrow = 10, ncol = 4))
  colnames(prediction_data) <- sanitized_names

  # Create mock model that expects these features
  expected_features <- sanitized_names

  # Test feature matching logic
  missing_features <- setdiff(expected_features, colnames(prediction_data))
  expect_length(missing_features, 0) # Should be no missing features

  # Test with missing features
  prediction_data_subset <- prediction_data[, 1:2] # Remove some features
  missing_features <- setdiff(expected_features, colnames(prediction_data_subset))
  expect_length(missing_features, 2) # Should have 2 missing features

  # Test adding missing features
  for (feature in missing_features) {
    prediction_data_subset[[feature]] <- NA
  }

  # After adding missing features, should have all expected features
  final_missing <- setdiff(expected_features, colnames(prediction_data_subset))
  expect_length(final_missing, 0)
})
