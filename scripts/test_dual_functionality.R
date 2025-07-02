#!/usr/bin/env Rscript

# Test Script for scCulturePredict Dual Functionality (Build/Predict Modes)
# This script tests the new dual functionality implemented in scumap()

# Load required libraries
library(scCulturePredict)
library(testthat)

# Set up test directories
test_base_dir <- "./test_dual_functionality"
training_dir <- file.path(test_base_dir, "training")
prediction_dir <- file.path(test_base_dir, "prediction")

# Clean up any existing test directories
if (dir.exists(test_base_dir)) {
  unlink(test_base_dir, recursive = TRUE)
}

# Create test directories
dir.create(test_base_dir, recursive = TRUE)
dir.create(training_dir, recursive = TRUE)
dir.create(prediction_dir, recursive = TRUE)

cat("=== Testing scCulturePredict Dual Functionality ===\n")
cat("Build Mode: Generate fingerprints from labeled data\n")
cat("Predict Mode: Apply fingerprints to unlabeled data\n\n")

# Test data paths (adjust these to your actual test data)
# Note: You'll need to provide actual paths to your test data
labeled_data_dir <- "./DATA_edit" # Replace with your labeled training data
unlabeled_data_dir <- "./DATA_edit" # Replace with your unlabeled test data
kegg_file <- "sce00001.keg" # Replace with your KEGG file path

# Check if test data exists
if (!dir.exists(labeled_data_dir)) {
  stop("Labeled training data directory not found: ", labeled_data_dir)
}
if (!dir.exists(unlabeled_data_dir)) {
  stop("Unlabeled test data directory not found: ", unlabeled_data_dir)
}
if (!file.exists(kegg_file)) {
  stop("KEGG file not found: ", kegg_file)
}

cat("Test data paths validated ✓\n\n")

# =============================================================================
# TEST 1: BUILD MODE - Generate Fingerprints
# =============================================================================

cat("=== TEST 1: BUILD MODE ===\n")
cat("Generating fingerprints from labeled training data...\n")

start_time <- Sys.time()

tryCatch(
  {
    training_results <- scumap(
      data_dir = labeled_data_dir,
      kegg_file = kegg_file,
      output_dir = training_dir,
      mode = "build",
      experiment_id = "test_training",
      progress = TRUE,
      verbose = TRUE
    )

    cat("Build mode completed successfully ✓\n")

    # Validate build mode results
    cat("Validating build mode results...\n")

    # Check that all expected components are returned
    expected_components <- c(
      "seurat_object", "pathway_results", "prediction_results",
      "evaluation_results", "fingerprint_file"
    )
    missing_components <- setdiff(expected_components, names(training_results))

    if (length(missing_components) > 0) {
      stop("Missing components in build mode results: ", paste(missing_components, collapse = ", "))
    }

    # Check that fingerprint file was created
    if (!file.exists(training_results$fingerprint_file)) {
      stop("Fingerprint file was not created: ", training_results$fingerprint_file)
    }

    # Check that Seurat object has predictions
    seurat_obj <- training_results$seurat_object
    if (is.null(seurat_obj$classification_pred)) {
      stop("Seurat object missing classification predictions")
    }

    # Check evaluation results
    eval_results <- training_results$evaluation_results
    if (is.null(eval_results$overall_accuracy)) {
      stop("Missing overall accuracy in evaluation results")
    }

    cat("Build mode validation passed ✓\n")
    cat("Fingerprint file created: ", training_results$fingerprint_file, "\n")
    cat("Training accuracy: ", round(eval_results$overall_accuracy, 3), "\n")
    cat("Number of cells trained: ", ncol(seurat_obj), "\n\n")
  },
  error = function(e) {
    cat("Build mode FAILED: ", e$message, "\n")
    stop("Build mode test failed")
  }
)

build_time <- Sys.time() - start_time
cat("Build mode time: ", round(as.numeric(build_time), 2), " seconds\n\n")

# =============================================================================
# TEST 2: PREDICT MODE - Apply Fingerprints
# =============================================================================

cat("=== TEST 2: PREDICT MODE ===\n")
cat("Applying fingerprints to unlabeled data...\n")

start_time <- Sys.time()

tryCatch(
  {
    prediction_results <- scumap(
      data_dir = unlabeled_data_dir,
      output_dir = prediction_dir,
      mode = "predict",
      fingerprint_file = training_results$fingerprint_file,
      experiment_id = "test_predictions",
      progress = TRUE,
      verbose = TRUE
    )

    cat("Predict mode completed successfully ✓\n")

    # Validate predict mode results
    cat("Validating predict mode results...\n")

    # Check that all expected components are returned
    expected_components <- c(
      "seurat_object", "pathway_results", "prediction_results",
      "evaluation_results", "fingerprint_source"
    )
    missing_components <- setdiff(expected_components, names(prediction_results))

    if (length(missing_components) > 0) {
      stop("Missing components in predict mode results: ", paste(missing_components, collapse = ", "))
    }

    # Check that Seurat object has predictions and confidence scores
    pred_seurat <- prediction_results$seurat_object
    if (is.null(pred_seurat$classification_pred)) {
      stop("Seurat object missing classification predictions")
    }
    if (is.null(pred_seurat$prediction_confidence)) {
      stop("Seurat object missing prediction confidence scores")
    }

    # Check that confidence scores are in valid range [0, 1]
    confidence_scores <- pred_seurat$prediction_confidence
    if (any(confidence_scores < 0 | confidence_scores > 1, na.rm = TRUE)) {
      stop("Confidence scores outside valid range [0, 1]")
    }

    cat("Predict mode validation passed ✓\n")
    cat("Number of predictions made: ", ncol(pred_seurat), "\n")
    cat("Mean confidence score: ", round(mean(confidence_scores, na.rm = TRUE), 3), "\n")
    cat("Prediction summary:\n")
    print(table(pred_seurat$classification_pred))
  },
  error = function(e) {
    cat("Predict mode FAILED: ", e$message, "\n")
    stop("Predict mode test failed")
  }
)

predict_time <- Sys.time() - start_time
cat("Predict mode time: ", round(as.numeric(predict_time), 2), " seconds\n\n")

# =============================================================================
# TEST 3: ERROR HANDLING - Invalid Parameters
# =============================================================================

cat("=== TEST 3: ERROR HANDLING ===\n")

# Test invalid mode
cat("Testing invalid mode parameter...\n")
expect_error(
  scumap(data_dir = labeled_data_dir, output_dir = training_dir, mode = "invalid"),
  "mode must be either 'build' or 'predict'"
)
cat("Invalid mode error handling ✓\n")

# Test build mode without KEGG file
cat("Testing build mode without KEGG file...\n")
expect_error(
  scumap(data_dir = labeled_data_dir, output_dir = training_dir, mode = "build"),
  "kegg_file must be provided"
)
cat("Build mode KEGG validation ✓\n")

# Test predict mode without fingerprint file
cat("Testing predict mode without fingerprint file...\n")
expect_error(
  scumap(data_dir = unlabeled_data_dir, output_dir = prediction_dir, mode = "predict"),
  "fingerprint_file must be provided"
)
cat("Predict mode fingerprint validation ✓\n")

# Test predict mode with invalid fingerprint file
cat("Testing predict mode with invalid fingerprint file...\n")
expect_error(
  scumap(
    data_dir = unlabeled_data_dir,
    output_dir = prediction_dir,
    mode = "predict",
    fingerprint_file = "nonexistent_file.rds"
  ),
  "Fingerprint file not found"
)
cat("Invalid fingerprint file error handling ✓\n")

cat("Error handling tests passed ✓\n\n")

# =============================================================================
# TEST 4: FILE OUTPUTS - Validate Generated Files
# =============================================================================

cat("=== TEST 4: FILE OUTPUTS ===\n")

# Check build mode outputs
cat("Checking build mode output files...\n")
build_files <- c(
  "seurat_object.rds",
  "pathway_results.rds",
  "prediction_results.rds",
  "evaluation_results.rds",
  "scCulturePredict_fingerprints.rds"
)

for (file in build_files) {
  file_path <- file.path(training_dir, file)
  if (!file.exists(file_path)) {
    stop("Missing build mode output file: ", file_path)
  }
}
cat("Build mode output files validated ✓\n")

# Check predict mode outputs
cat("Checking predict mode output files...\n")
predict_files <- c(
  "predicted_seurat_object.rds",
  "applied_pathway_results.rds",
  "prediction_results.rds",
  "prediction_evaluation.rds"
)

for (file in predict_files) {
  file_path <- file.path(prediction_dir, file)
  if (!file.exists(file_path)) {
    stop("Missing predict mode output file: ", file_path)
  }
}
cat("Predict mode output files validated ✓\n")

# Test loading saved files
cat("Testing loading of saved files...\n")
tryCatch(
  {
    # Load build mode results
    saved_seurat <- readRDS(file.path(training_dir, "seurat_object.rds"))
    saved_fingerprints <- readRDS(file.path(training_dir, "scCulturePredict_fingerprints.rds"))

    # Load predict mode results
    saved_predictions <- readRDS(file.path(prediction_dir, "predicted_seurat_object.rds"))

    cat("File loading successful ✓\n")
  },
  error = function(e) {
    stop("Failed to load saved files: ", e$message)
  }
)

# =============================================================================
# TEST 5: FINGERPRINT COMPATIBILITY
# =============================================================================

cat("=== TEST 5: FINGERPRINT COMPATIBILITY ===\n")

# Load and validate fingerprint file structure
cat("Validating fingerprint file structure...\n")
fingerprint_data <- readRDS(training_results$fingerprint_file)

required_components <- c("kegg_pathways", "pathway_results", "svm_model", "similarity_model", "metadata")
missing_components <- setdiff(required_components, names(fingerprint_data))

if (length(missing_components) > 0) {
  stop("Invalid fingerprint structure. Missing: ", paste(missing_components, collapse = ", "))
}

# Check metadata
metadata <- fingerprint_data$metadata
if (is.null(metadata$created_date) || is.null(metadata$training_accuracy)) {
  stop("Missing required metadata in fingerprint file")
}

cat("Fingerprint file structure validated ✓\n")
cat("Metadata - Created: ", as.character(metadata$created_date), "\n")
cat("Metadata - Training accuracy: ", round(metadata$training_accuracy, 3), "\n")
cat("Metadata - Number of pathways: ", metadata$n_pathways, "\n")

# =============================================================================
# TEST SUMMARY
# =============================================================================

cat("\n=== TEST SUMMARY ===\n")
cat("✓ Build mode functionality\n")
cat("✓ Predict mode functionality\n")
cat("✓ Error handling and validation\n")
cat("✓ File output generation\n")
cat("✓ Fingerprint compatibility\n")
cat("✓ All tests PASSED\n\n")

total_time <- build_time + predict_time
cat("Total test time: ", round(as.numeric(total_time), 2), " seconds\n")

# =============================================================================
# CLEANUP (Optional)
# =============================================================================

cat("\n=== CLEANUP ===\n")
cat("Test files are saved in: ", test_base_dir, "\n")
cat("To clean up test files, run: unlink('", test_base_dir, "', recursive = TRUE)\n")

# Uncomment the following line to automatically clean up test files
# unlink(test_base_dir, recursive = TRUE)
# cat("Test files cleaned up ✓\n")

cat("\n=== DUAL FUNCTIONALITY TESTING COMPLETED ===\n")
cat("The scCulturePredict package dual functionality (build/predict modes) is working correctly!\n")
