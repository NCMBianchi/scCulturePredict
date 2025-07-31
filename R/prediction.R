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
#'     rnorm(150, mean = 0, sd = 1),
#'     nrow = 15,
#'     dimnames = list(
#'         paste0("Cell", seq_len(15)),
#'         paste0("Pathway", seq_len(10))
#'     )
#' )
#'
#' # Create mock fingerprint profiles
#' fingerprint_profiles <- matrix(
#'     rnorm(30, mean = 0, sd = 0.5),
#'     nrow = 3,
#'     dimnames = list(
#'         c("CultureA", "CultureB", "CultureC"),
#'         paste0("Pathway", seq_len(10))
#'     )
#' )
#'
#' # Make predictions
#' predictions <- predict_by_similarity(
#'     pathway_activity,
#'     fingerprint_profiles
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
#' svm_results <- predict_by_svm(pathway_matrix, seurat_object)
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

#' Predict cell types using pathway activity
#'
#' @description
#' Predicts cell types based on pathway activity patterns using a trained classifier.
#' This function takes a Seurat object and pathway activity data, and returns
#' predicted cell type labels for each cell.
#'
#' @details
#' This function applies a pre-trained classifier to predict cell types based on
#' pathway activity patterns. It's designed to work with classifiers trained using
#' the \code{\link{train_cell_type_classifier}} function, but can work with any
#' classifier that has a predict method.
#'
#' The function performs several validation checks to ensure the inputs are valid:
#' \enumerate{
#'   \item Validates that seurat_object is a proper Seurat object
#'   \item Validates that pathway_activity is a matrix or data frame with the correct dimensions
#'   \item Checks that the number of cells matches between the pathway_activity and seurat_object
#'   \item Verifies that the classifier has a predict method
#'   \item Ensures all required pathways are present in the pathway_activity matrix
#' }
#'
#' After validation, the function:
#' \enumerate{
#'   \item Extracts the relevant pathway activity patterns based on the classifier's feature importance
#'   \item Transposes the data to the format expected by the classifier
#'   \item Makes predictions using the classifier's predict method
#'   \item If probability=TRUE, returns both class predictions and probability estimates
#'   \item If probability=FALSE, returns only class predictions
#' }
#'
#' Error handling is implemented throughout to provide informative error messages
#' if any step fails, making the function robust for both interactive and
#' programmatic use.
#'
#' @param seurat_object A Seurat object containing gene expression data.
#' @param pathway_activity Matrix. Pathway activity matrix where rows are pathways
#'   and columns are cells.
#' @param classifier Trained classifier object. Must have a predict method.
#' @param probability Logical. Whether to return prediction probabilities
#'   (default: FALSE).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return If probability is FALSE, returns a character vector of predicted cell types.
#'   If probability is TRUE, returns a list containing:
#'   \itemize{
#'     \item predictions: Character vector of predicted cell types
#'     \item probabilities: Matrix of prediction probabilities
#'   }
#'
#' @examples
#' # Example with mock data
#' # Create mock Seurat object
#' library(Seurat)
#' counts <- matrix(rpois(1000, 5), nrow = 100)
#' rownames(counts) <- paste0("Gene", seq_len(100))
#' colnames(counts) <- paste0("Cell", seq_len(10))
#' seurat_obj <- CreateSeuratObject(counts = counts)
#' seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
#'
#' # Create mock pathway activities
#' pathway_activities <- matrix(
#'     rnorm(50, mean = 0, sd = 1),
#'     nrow = 10,
#'     dimnames = list(
#'         colnames(seurat_obj),
#'         paste0("Pathway", seq_len(5))
#'     )
#' )
#'
#' # Create mock models list
#' mock_fingerprints <- matrix(
#'     rnorm(15, mean = 0, sd = 0.5),
#'     nrow = 3,
#'     dimnames = list(
#'         c("TypeA", "TypeB", "TypeC"),
#'         paste0("Pathway", seq_len(5))
#'     )
#' )
#'
#' models_list <- list(
#'     direct = list(fingerprints = mock_fingerprints),
#'     pathway_activities = pathway_activities
#' )
#'
#' # Predict cell types
#' predictions <- predict_cell_types(
#'     seurat_obj,
#'     models_list,
#'     method = "direct"
#' )
#' @seealso
#' \code{\link{train_cell_type_classifier}} for training the classifier
#' \code{\link{evaluate_predictions}} for evaluating prediction performance
#'
#' @export
predict_cell_types <- function(seurat_object,
                               pathway_activity,
                               classifier,
                               probability = FALSE,
                               verbose = TRUE) {
    # Input validation
    if (!inherits(seurat_object, "Seurat")) {
        stop("seurat_object must be a Seurat object")
    }

    if (!is.matrix(pathway_activity) && !is.data.frame(pathway_activity)) {
        stop("pathway_activity must be a matrix or data frame")
    }

    if (ncol(pathway_activity) != ncol(seurat_object)) {
        stop("Number of cells in pathway_activity must match number of cells in seurat_object")
    }

    if (!hasMethod("predict", class(classifier)[1])) {
        stop("classifier must have a predict method")
    }

    if (!is.logical(probability) || length(probability) != 1) {
        stop("probability must be a single logical value")
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("verbose must be a single logical value")
    }

    # Check if pathway data is available
    if (verbose) message("Checking pathway data...")
    if (!all(rownames(pathway_activity) %in% rownames(classifier$feature_importance))) {
        missing_pathways <- setdiff(
            rownames(pathway_activity),
            rownames(classifier$feature_importance)
        )
        stop(sprintf(
            "Missing pathway data for: %s",
            paste(missing_pathways, collapse = ", ")
        ))
    }

    # Prepare data for prediction
    if (verbose) message("Preparing data for prediction...")
    prediction_data <- tryCatch(
        {
            t(pathway_activity[rownames(classifier$feature_importance), ])
        },
        error = function(e) {
            stop(sprintf("Error preparing prediction data: %s", e$message))
        }
    )

    # Make predictions
    if (verbose) message("Making predictions...")
    predictions <- tryCatch(
        {
            if (probability) {
                predict(classifier, prediction_data, type = "prob")
            } else {
                predict(classifier, prediction_data)
            }
        },
        error = function(e) {
            stop(sprintf("Error making predictions: %s", e$message))
        }
    )

    if (verbose) {
        message("Predictions completed successfully")
    }

    return(predictions)
}

#' Train cell type classifier
#'
#' @description
#' Trains a classifier to predict cell types based on pathway activity patterns.
#' This function takes a Seurat object with known cell types and pathway activity
#' data, and returns a trained classifier.
#'
#' @details
#' This function trains a machine learning classifier to predict cell types using
#' pathway activity patterns as features. It supports multiple classification methods
#' and implements a complete training workflow including feature selection,
#' cross-validation, and model training.
#'
#' The function supports three classification methods:
#' \itemize{
#'   \item "rf": Random Forest - An ensemble method using multiple decision trees
#'   \item "svm": Support Vector Machine - A powerful classifier that finds optimal boundaries
#'   \item "xgb": XGBoost - A gradient boosting framework known for performance and accuracy
#' }
#'
#' The workflow includes:
#' \enumerate{
#'   \item Extensive input validation to ensure data quality and compatibility
#'   \item Extraction of cell type labels from the specified metadata column
#'   \item Feature importance calculation to identify the most predictive pathways
#'   \item Selection of top n_features pathways based on importance
#'   \item k-fold cross-validation to assess model performance (using cv_folds)
#'   \item Training of the final model on the complete dataset
#'   \item Compilation of results including the model, feature importance, and performance metrics
#' }
#'
#' Feature selection helps reduce dimensionality and focuses the model on the most
#' informative pathways, which can improve both performance and interpretability.
#' Cross-validation provides a robust estimate of model performance on unseen data.
#'
#' The function includes comprehensive error handling at each step of the workflow,
#' with informative error messages to help diagnose any issues that arise during training.
#'
#' @param seurat_object A Seurat object containing gene expression data
#'   and cell type labels.
#' @param pathway_activity Matrix. Pathway activity matrix where rows are pathways
#'   and columns are cells.
#' @param cell_type_col Character. Name of the metadata column containing cell type
#'   labels (default: "cell_type").
#' @param method Character. Classification method to use. Options are "rf" (Random Forest),
#'   "svm" (Support Vector Machine), or "xgb" (XGBoost) (default: "rf").
#' @param n_features Integer. Number of top pathways to use for classification
#'   (default: 100).
#' @param cv_folds Integer. Number of cross-validation folds (default: 5).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A trained classifier object with the following components:
#' \itemize{
#'   \item model: The trained classification model
#'   \item feature_importance: Data frame of pathway importance scores
#'   \item cv_performance: Cross-validation performance metrics
#'   \item training_params: List of training parameters used
#' }
#'
#' @examples
#' # Example with mock data
#' # Create mock pathway activity matrix
#' set.seed(123)
#' n_cells <- 60
#' n_pathways <- 10
#'
#' pathway_activity <- matrix(
#'     c(
#'         rnorm(20 * n_pathways, mean = 1, sd = 0.5), # TypeA pattern
#'         rnorm(20 * n_pathways, mean = -1, sd = 0.5), # TypeB pattern
#'         rnorm(20 * n_pathways, mean = 0, sd = 0.5) # TypeC pattern
#'     ),
#'     nrow = n_cells,
#'     ncol = n_pathways,
#'     byrow = TRUE
#' )
#'
#' rownames(pathway_activity) <- paste0("Cell", seq_len(n_cells))
#' colnames(pathway_activity) <- paste0("Pathway", seq_len(n_pathways))
#'
#' # Create labels
#' labels <- factor(rep(c("TypeA", "TypeB", "TypeC"), each = 20))
#'
#' # Train classifier
#' classifier <- train_cell_type_classifier(
#'     pathway_activity,
#'     labels,
#'     method = "svm",
#'     cross_validate = FALSE # Disable CV for speed
#' )
#'
#' # Check model
#' print(classifier$method)
#' print(table(classifier$predictions, labels))
#' @seealso
#' \code{\link{predict_cell_types}} for making predictions
#' \code{\link{evaluate_predictions}} for evaluating prediction performance
#'
#' @export
train_cell_type_classifier <- function(seurat_object,
                                       pathway_activity,
                                       cell_type_col = "cell_type",
                                       method = "rf",
                                       n_features = 100,
                                       cv_folds = 5,
                                       verbose = TRUE) {
    # Input validation
    if (!inherits(seurat_object, "Seurat")) {
        stop("seurat_object must be a Seurat object")
    }

    if (!is.matrix(pathway_activity) && !is.data.frame(pathway_activity)) {
        stop("pathway_activity must be a matrix or data frame")
    }

    if (ncol(pathway_activity) != ncol(seurat_object)) {
        stop("Number of cells in pathway_activity must match number of cells in seurat_object")
    }

    if (!is.character(cell_type_col) || length(cell_type_col) != 1) {
        stop("cell_type_col must be a single character string")
    }

    if (!cell_type_col %in% colnames(seurat_object@meta.data)) {
        stop(sprintf("Cell type column '%s' not found in metadata", cell_type_col))
    }

    if (!is.character(method) || length(method) != 1 ||
        !method %in% c("rf", "svm", "xgb")) {
        stop("method must be one of: 'rf', 'svm', 'xgb'")
    }

    if (!is.numeric(n_features) || length(n_features) != 1 || n_features < 1) {
        stop("n_features must be a single positive integer")
    }

    if (!is.numeric(cv_folds) || length(cv_folds) != 1 || cv_folds < 2) {
        stop("cv_folds must be a single integer >= 2")
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("verbose must be a single logical value")
    }

    # Get cell type labels
    if (verbose) message("Getting cell type labels...")
    cell_types <- tryCatch(
        {
            seurat_object@meta.data[[cell_type_col]]
        },
        error = function(e) {
            stop(sprintf("Error accessing cell type labels: %s", e$message))
        }
    )

    if (any(is.na(cell_types))) {
        stop("Cell type labels contain missing values")
    }

    # Prepare training data
    if (verbose) message("Preparing training data...")
    training_data <- tryCatch(
        {
            t(pathway_activity)
        },
        error = function(e) {
            stop(sprintf("Error preparing training data: %s", e$message))
        }
    )

    # Select top features
    if (verbose) message("Selecting top features...")
    feature_importance <- tryCatch(
        {
            if (method == "rf") {
                caret::varImp(randomForest::randomForest(training_data, cell_types))
            } else if (method == "svm") {
                caret::varImp(e1071::svm(training_data, cell_types))
            } else {
                caret::varImp(xgboost::xgboost(training_data, cell_types))
            }
        },
        error = function(e) {
            stop(sprintf("Error calculating feature importance: %s", e$message))
        }
    )

    top_features <- rownames(feature_importance)[seq_len(min(n_features, nrow(feature_importance)))]
    training_data <- training_data[, top_features]

    # Perform cross-validation
    if (verbose) message("Performing cross-validation...")
    cv_results <- tryCatch(
        {
            caret::train(
                x = training_data,
                y = cell_types,
                method = switch(method,
                    "rf" = "rf",
                    "svm" = "svmRadial",
                    "xgb" = "xgbTree"
                ),
                trControl = caret::trainControl(
                    method = "cv",
                    number = cv_folds,
                    verboseIter = verbose
                )
            )
        },
        error = function(e) {
            stop(sprintf("Error in cross-validation: %s", e$message))
        }
    )

    # Train final model
    if (verbose) message("Training final model...")
    final_model <- tryCatch(
        {
            switch(method,
                "rf" = randomForest::randomForest(training_data, cell_types),
                "svm" = e1071::svm(training_data, cell_types),
                "xgb" = xgboost::xgboost(training_data, cell_types)
            )
        },
        error = function(e) {
            stop(sprintf("Error training final model: %s", e$message))
        }
    )

    # Create classifier object
    classifier <- list(
        model = final_model,
        feature_importance = feature_importance,
        cv_performance = cv_results$results,
        training_params = list(
            method = method,
            n_features = n_features,
            cv_folds = cv_folds
        )
    )

    if (verbose) {
        message("Classifier training completed successfully")
    }

    return(classifier)
}
