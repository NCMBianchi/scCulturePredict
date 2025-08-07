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
  # Normalize expression matrix
  pathway_matrix_norm <- scale(pathway_matrix)

  # Calculate cosine similarity
  cosine_similarity <- function(v1, v2) {
    valid_idx <- !is.na(v1) & !is.na(v2)
    sum(v1[valid_idx] * v2[valid_idx]) /
      (sqrt(sum(v1[valid_idx]^2)) * sqrt(sum(v2[valid_idx]^2)))
  }

  # Build similarity matrix
  similarity_matrix <- matrix(nrow = nrow(pathway_matrix_norm), ncol = ncol(signature_matrix))
  for (i in seq_len(nrow(pathway_matrix_norm))) {
    for (j in seq_len(ncol(signature_matrix))) {
      similarity_matrix[i, j] <- cosine_similarity(pathway_matrix_norm[i, ], signature_matrix[, j])
    }
  }
  colnames(similarity_matrix) <- colnames(signature_matrix)

  # Direct prediction (highest similarity)
  predicted_direct <- colnames(signature_matrix)[apply(similarity_matrix, 1, which.max)]

  # Threshold prediction (only if above threshold)
  predicted_threshold <- apply(similarity_matrix, 1, function(row) {
    if (max(row) > threshold) colnames(signature_matrix)[which.max(row)] else NA
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
  # Prepare data for classification
  pathway_df <- as.data.frame(pathway_matrix[, colSums(is.na(pathway_matrix)) == 0])

  # Sanitize column names to ensure they are valid R variable names
  # This is critical for SVM model consistency between training and prediction
  colnames(pathway_df) <- make.names(colnames(pathway_df))

  combined_data <- cbind(pathway_df, Condition = seurat_object@meta.data$sample)

  # Train-test split
  train_idx <- sample(seq_len(nrow(combined_data)), size = 0.8 * nrow(combined_data))

  # Scale data
  training_data <- combined_data[train_idx, ]
  scaled_train <- scale(training_data[, -ncol(training_data)])
  scaled_train <- as.data.frame(scaled_train)
  scaled_train$Condition <- as.factor(training_data$Condition)

  testing_data <- combined_data[-train_idx, ]
  scaled_test <- scale(testing_data[, -ncol(testing_data)])
  scaled_test <- as.data.frame(scaled_test)
  scaled_test$Condition <- as.factor(testing_data$Condition)

  # Train SVM model
  svm_model <- e1071::svm(Condition ~ ., data = scaled_train, kernel = "radial")

  # Make predictions on full dataset
  scaled_all <- scale(pathway_df)
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
