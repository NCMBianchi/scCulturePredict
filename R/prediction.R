#' Predict growth media by similarity
#'
#' @description
#' Predicts growth media conditions based on cosine similarity between pathway expression
#' and condition-specific signatures.
#'
#' @details
#' This function implements a similarity-based approach for predicting growth media conditions
#' by comparing pathway expression patterns against known signatures. It uses cosine similarity
#' as the metric for comparison, which measures the cosine of the angle between two vectors.
#'
#' The function performs the following steps:
#' \enumerate{
#'   \item Normalizes the input pathway expression matrix to standardize values
#'   \item Calculates the cosine similarity between each cell's pathway expression and each condition signature
#'   \item Creates a similarity matrix with cells as rows and conditions as columns
#'   \item Makes direct predictions by selecting the condition with highest similarity for each cell
#'   \item Makes threshold-based predictions, returning NA if maximum similarity is below threshold
#' }
#'
#' Cosine similarity ranges from -1 (completely opposite) to 1 (identical direction),
#' with 0 indicating orthogonality (no relationship). The threshold parameter allows
#' control over prediction stringency - higher values require stronger similarity for prediction.
#'
#' @param pathway_matrix Matrix of pathway expression values.
#' @param signature_matrix Matrix of condition-specific signature profiles.
#' @param threshold Numeric threshold for similarity-based predictions. Default is 0.1.
#' @param verbose Logical indicating whether to print progress messages. Default is TRUE.
#'
#' @return A list containing:
#' \itemize{
#'   \item similarity_matrix: Matrix of cosine similarities
#'   \item predicted_direct: Direct predictions (highest similarity)
#'   \item predicted_threshold: Threshold-based predictions
#' }
#' @export
#'
#' @examples
#' # Example with mock data
#' # Create mock pathway activity matrix
#' pathway_activity <- matrix(
#'   rnorm(150, mean = 0, sd = 1),
#'   nrow = 15,
#'   dimnames = list(
#'     paste0("Cell", seq_len(15)),
#'     paste0("Pathway", seq_len(10))
#'   )
#' )
#'
#' # Create mock fingerprint profiles
#' # Note: fingerprints should have pathways as rows and cultures as columns
#' fingerprint_profiles <- matrix(
#'   rnorm(30, mean = 0, sd = 0.5),
#'   nrow = 10,
#'   dimnames = list(
#'     paste0("Pathway", seq_len(10)),
#'     c("CultureA", "CultureB", "CultureC")
#'   )
#' )
#'
#' # Make predictions
#' predictions <- predict_by_similarity(
#'   pathway_activity,
#'   fingerprint_profiles
#' )
#'
#' # Check results
#' print(head(predictions))
predict_by_similarity <- function(pathway_matrix, signature_matrix, threshold = 0.1, verbose = TRUE) {
  # Helper function for safe scaling that handles zero variance
  safe_scale <- function(x) {
    scaled <- scale(x)
    # Replace NaN values (from zero variance columns) with 0
    scaled[is.nan(scaled)] <- 0
    return(scaled)
  }

  # Find common pathways between the two matrices
  pathway_names_new <- colnames(pathway_matrix)
  pathway_names_ref <- rownames(signature_matrix)
  common_pathways <- intersect(pathway_names_new, pathway_names_ref)

  if (length(common_pathways) == 0) {
    warning("No common pathways found between new data and reference fingerprint")
    # Return empty predictions
    return(list(
      similarity_matrix = matrix(NA, nrow = nrow(pathway_matrix), ncol = ncol(signature_matrix)),
      predicted_direct = rep(NA, nrow(pathway_matrix)),
      predicted_threshold = rep(NA, nrow(pathway_matrix))
    ))
  }

  if (verbose) {
    message(sprintf(
      "Using %d common pathways out of %d in reference and %d in new data",
      length(common_pathways),
      length(pathway_names_ref),
      length(pathway_names_new)
    ))
  }

  # Subset both matrices to use only common pathways
  pathway_matrix_subset <- pathway_matrix[, common_pathways, drop = FALSE]
  signature_matrix_subset <- signature_matrix[common_pathways, , drop = FALSE]

  # Normalize expression matrix (using subset)
  pathway_matrix_norm <- safe_scale(pathway_matrix_subset)

  # Calculate cosine similarity
  cosine_similarity <- function(v1, v2) {
    valid_idx <- !is.na(v1) & !is.na(v2) & !is.nan(v1) & !is.nan(v2)
    if (sum(valid_idx) == 0) {
      return(0)
    } # No valid values

    v1_valid <- v1[valid_idx]
    v2_valid <- v2[valid_idx]

    norm1 <- sqrt(sum(v1_valid^2))
    norm2 <- sqrt(sum(v2_valid^2))

    # Handle zero vectors
    if (norm1 == 0 || norm2 == 0) {
      return(0)
    }

    return(sum(v1_valid * v2_valid) / (norm1 * norm2))
  }

  # Build similarity matrix (using subsets)
  similarity_matrix <- matrix(nrow = nrow(pathway_matrix_norm), ncol = ncol(signature_matrix_subset))
  for (i in seq_len(nrow(pathway_matrix_norm))) {
    for (j in seq_len(ncol(signature_matrix_subset))) {
      similarity_matrix[i, j] <- cosine_similarity(pathway_matrix_norm[i, ], signature_matrix_subset[, j])
    }
  }
  colnames(similarity_matrix) <- colnames(signature_matrix_subset)

  # Direct prediction (highest similarity)
  predicted_direct <- colnames(signature_matrix_subset)[apply(similarity_matrix, 1, which.max)]

  # Threshold prediction (only if above threshold)
  predicted_threshold <- apply(similarity_matrix, 1, function(row) {
    valid_row <- row[!is.na(row) & !is.nan(row)]
    if (length(valid_row) == 0) {
      return(NA)
    }
    max_val <- max(valid_row, na.rm = TRUE)
    if (is.na(max_val) || is.nan(max_val) || max_val <= threshold) {
      return(NA)
    }
    return(colnames(signature_matrix_subset)[which.max(row)])
  })

  return(list(
    similarity_matrix = similarity_matrix,
    predicted_direct = predicted_direct,
    predicted_threshold = predicted_threshold
  ))
}

#' Predict growth media using SVM
#'
#' @description
#' Predicts growth media conditions using Support Vector Machine (SVM) classification.
#'
#' @details
#' This function implements a machine learning approach for predicting growth media
#' conditions using Support Vector Machine (SVM) classification. Unlike similarity-based
#' methods, SVM can capture more complex, non-linear relationships between pathway
#' expression patterns and growth conditions.
#'
#' The workflow consists of:
#' \enumerate{
#'   \item Preparing data by removing NA values and combining with condition labels
#'   \item Splitting data into training (80%) and testing (20%) sets
#'   \item Scaling features to standardize the range and distribution
#'   \item Training an SVM model with radial basis function kernel
#'   \item Making predictions on the full dataset
#'   \item Evaluating model performance on the test set
#' }
#'
#' The radial basis function kernel allows the SVM to capture non-linear relationships
#' in the data, potentially improving prediction accuracy over linear methods. The
#' function handles data preparation, model training, prediction, and evaluation
#' in an integrated workflow.
#'
#' @param pathway_matrix Matrix of pathway expression values.
#' @param seurat_object A Seurat object containing single-cell data.
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#' \itemize{
#'   \item predictions: SVM predictions for all cells
#'   \item svm_model: Trained SVM model
#'   \item confusion_matrix: Confusion matrix for test set
#'   \item accuracy: Overall accuracy on test set
#' }
#' @export
#'
#' @examples
#' # Create mock data with sufficient samples for SVM
#' library(Seurat)
#' set.seed(123)
#'
#' # Create a mock Seurat object with sample metadata
#' # Need at least 20 cells for proper train/test split
#' counts <- matrix(rpois(1000, 5), nrow = 50, ncol = 20)
#' rownames(counts) <- paste0("Gene", seq_len(50))
#' colnames(counts) <- paste0("Cell", seq_len(20))
#'
#' # Create metadata with proper row names
#' metadata <- data.frame(
#'   row.names = colnames(counts),
#'   sample = rep(c("SampleA", "SampleB"), each = 10)
#' )
#' seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
#'
#' # Create a mock pathway activity matrix directly
#' # This represents pathway activities for each cell
#' pathway_matrix <- matrix(
#'   rnorm(20 * 5, mean = 0, sd = 1), # 20 cells x 5 pathways
#'   nrow = 20,
#'   ncol = 5,
#'   dimnames = list(
#'     colnames(counts), # Cell names must match Seurat object
#'     paste0("Pathway", seq_len(5)) # Pathway names
#'   )
#' )
#'
#' # Make pathway activities slightly different between samples
#' # to ensure SVM can find patterns
#' pathway_matrix[1:10, ] <- pathway_matrix[1:10, ] + 0.5 # SampleA cells
#' pathway_matrix[11:20, ] <- pathway_matrix[11:20, ] - 0.5 # SampleB cells
#'
#' # Run SVM prediction
#' svm_results <- predict_by_svm(
#'   pathway_matrix,
#'   seurat_obj,
#'   verbose = FALSE
#' )
#'
#' # Check results
#' print(svm_results$accuracy)
predict_by_svm <- function(pathway_matrix, seurat_object, verbose = TRUE) {
  # Filter out pathways with all NA values
  non_na_cols <- colSums(is.na(pathway_matrix)) < nrow(pathway_matrix)
  pathway_df <- as.data.frame(pathway_matrix[, non_na_cols, drop = FALSE])

  if (ncol(pathway_df) == 0) {
    stop("No valid pathways found in the data (all pathways have only NA values)")
  }

  # Store original pathway names before sanitizing
  original_pathway_names <- colnames(pathway_df)

  # Sanitize column names to ensure they are valid R variable names
  # This is critical for SVM model consistency between training and prediction
  colnames(pathway_df) <- make.names(colnames(pathway_df))

  # Create mapping between original and sanitized names for later use
  pathway_name_mapping <- setNames(colnames(pathway_df), original_pathway_names)

  combined_data <- cbind(pathway_df, Condition = seurat_object@meta.data$sample)

  # Validate data before train-test split
  n_samples <- nrow(combined_data)
  unique_conditions <- unique(combined_data$Condition)
  n_conditions <- length(unique_conditions)

  if (n_samples < 10) {
    stop("Insufficient data for SVM training: need at least 10 cells, got ", n_samples)
  }

  if (n_conditions < 2) {
    stop("SVM requires at least 2 different conditions, got ", n_conditions)
  }

  # Check minimum samples per condition
  min_samples_per_condition <- min(table(combined_data$Condition))
  if (min_samples_per_condition < 2) {
    stop("Each condition must have at least 2 cells for SVM training")
  }

  if (verbose && ncol(pathway_df) < 10) {
    message(sprintf("Note: Training SVM with only %d pathways (low pathway coverage)", ncol(pathway_df)))
  }

  # Train-test split with validation
  train_size <- floor(0.8 * n_samples)
  if (train_size < n_conditions * 2) {
    # If training set would be too small, use leave-one-out or adjust split
    train_size <- max(n_conditions * 2, floor(0.9 * n_samples))
  }

  train_idx <- sample(seq_len(n_samples), size = train_size)

  # Helper function for safe scaling that handles zero variance
  safe_scale <- function(x) {
    scaled <- scale(x)
    # Replace NaN values (from zero variance columns) with 0
    scaled[is.nan(scaled)] <- 0
    return(scaled)
  }

  # Scale data
  training_data <- combined_data[train_idx, ]
  scaled_train <- safe_scale(training_data[, -ncol(training_data)])
  scaled_train <- as.data.frame(scaled_train)
  scaled_train$Condition <- as.factor(training_data$Condition)

  testing_data <- combined_data[-train_idx, ]
  scaled_test <- safe_scale(testing_data[, -ncol(testing_data)])
  scaled_test <- as.data.frame(scaled_test)
  scaled_test$Condition <- as.factor(testing_data$Condition)

  # Validate training data has all conditions
  if (length(unique(scaled_train$Condition)) != n_conditions) {
    # Re-sample to ensure all conditions are in training set
    for (attempt in 1:10) {
      train_idx <- sample(seq_len(n_samples), size = train_size)
      training_data <- combined_data[train_idx, ]
      if (length(unique(training_data$Condition)) == n_conditions) break
    }
    # Re-scale with new split
    scaled_train <- safe_scale(training_data[, -ncol(training_data)])
    scaled_train <- as.data.frame(scaled_train)
    scaled_train$Condition <- as.factor(training_data$Condition)
  }

  # Train SVM model
  svm_model <- e1071::svm(Condition ~ ., data = scaled_train, kernel = "radial")

  # Store the features used in training for future predictions
  svm_model$training_features <- colnames(pathway_df)
  svm_model$original_pathway_names <- original_pathway_names

  # Make predictions on full dataset
  scaled_all <- safe_scale(pathway_df)
  scaled_all <- as.data.frame(scaled_all)
  # Ensure column names are preserved after scaling
  colnames(scaled_all) <- colnames(pathway_df)
  predictions <- predict(svm_model, newdata = scaled_all)

  # Test accuracy
  # Ensure column names match training data
  colnames(scaled_test) <- c(colnames(pathway_df)[seq_len(ncol(scaled_test) - 1)], "Condition")
  test_predictions <- predict(svm_model, newdata = scaled_test)
  confusion_matrix <- table(Predicted = test_predictions, Actual = scaled_test$Condition)
  accuracy <- sum(diag(confusion_matrix)) / sum(confusion_matrix)

  return(list(
    predictions = predictions,
    svm_model = svm_model,
    confusion_matrix = confusion_matrix,
    accuracy = accuracy
  ))
}
