#' Single-Cell Culture Media Prediction Pipeline
#'
#' @description
#' Runs the complete scCulturePredict analysis pipeline with two modes: "build" to generate
#' fingerprints from labeled training data, or "predict" to apply existing fingerprints
#' to new unlabeled data for cell culture media prediction.
#'
#' @param data_dir Character string specifying the directory containing the data files.
#' @param kegg_file Character string specifying the path to the KEGG pathway file.
#'   Required for "build" mode, optional for "predict" mode if fingerprint_file contains KEGG data.
#' @param output_dir Character string specifying the directory for output files.
#' @param mode Character string specifying the analysis mode. Either "build" (generate
#'   fingerprints from labeled data) or "predict" (apply existing fingerprints to unlabeled data).
#'   Default is "build".
#' @param fingerprint_file Character string specifying the path to saved fingerprint/model file.
#'   Required for "predict" mode, ignored in "build" mode.
#' @param experiment_id Character string specifying the experiment ID prefix in filenames.
#'   Default is "experiment".
#' @param perform_tsne Logical indicating whether to perform t-SNE. Default is TRUE.
#' @param progress Logical indicating whether to show progress bar. Default is FALSE.
#' @param parallel Logical indicating whether to use parallel processing. Default is FALSE.
#' @param n_cores Integer specifying the number of cores to use for parallel processing.
#'   Default is NULL (uses all available cores).
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#' \itemize{
#'   \item seurat_object: The processed Seurat object with UMAP coordinates and predictions
#'   \item pathway_results: Results from KEGG pathway analysis (build mode) or fingerprint data (predict mode)
#'   \item prediction_results: Results from similarity and SVM prediction analysis
#'   \item evaluation_results: Results from prediction evaluation (build mode) or prediction confidence (predict mode)
#'   \item fingerprint_file: Path to saved fingerprint file (build mode only)
#' }
#'
#' @details
#' The scCulture function operates in two distinct modes:
#'
#' **Build Mode** (\code{mode = "build"}):
#' Generates transferable fingerprints from labeled training data:
#' \enumerate{
#'   \item Validates input parameters and file existence
#'   \item Creates output directory if it doesn't exist
#'   \item Loads and preprocesses labeled single-cell data
#'   \item Performs dimensionality reduction (PCA, UMAP, and optionally t-SNE)
#'   \item Parses KEGG pathways and builds transcriptomic fingerprints
#'   \item Trains similarity-based and SVM models
#'   \item Evaluates model performance with cross-validation
#'   \item Saves fingerprints and models for future predictions
#'   \item Creates publication-ready visualizations and saves results
#' }
#'
#' **Predict Mode** (\code{mode = "predict"}):
#' Applies existing fingerprints to predict unlabeled data:
#' \enumerate{
#'   \item Validates input parameters and loads fingerprint file
#'   \item Loads and preprocesses unlabeled single-cell data
#'   \item Performs dimensionality reduction consistent with training data
#'   \item Applies pre-built fingerprints to new data
#'   \item Makes predictions using trained similarity and SVM models
#'   \item Calculates prediction confidence scores
#'   \item Creates prediction visualizations
#' }
#'
#' The pipeline integrates several key analytical components:
#' \itemize{
#'   \item Data loading and preprocessing using Seurat framework
#'   \item Dimensionality reduction using PCA, UMAP, and optional t-SNE
#'   \item KEGG pathway analysis for biological interpretation
#'   \item Dual prediction approaches: similarity-based and machine learning (SVM)
#'   \item Comprehensive evaluation with statistical testing
#'   \item Professional visualization suite with customizable plots
#' }
#'
#' Progress tracking options:
#' \itemize{
#'   \item \code{progress = TRUE}: Shows detailed progress bar with completion percentages
#'   \item \code{verbose = TRUE}: Displays step-by-step progress messages
#'   \item \code{parallel = TRUE}: Enables parallel processing for intensive computations
#' }
#'
#' Output files automatically saved to \code{output_dir}:
#' \itemize{
#'   \item \code{seurat_object.rds}: Complete Seurat object with all analysis results
#'   \item \code{pathway_results.rds}: KEGG pathway analysis results
#'   \item \code{prediction_results.rds}: Similarity and SVM prediction results
#'   \item \code{evaluation_results.rds}: Performance evaluation metrics
#'   \item \code{umap_plots.pdf}: UMAP visualization plots
#'   \item \code{evaluation_plots.pdf}: Prediction performance plots
#' }
#'
#' @examples
#' \dontrun{
#' # Build mode - generate fingerprints from labeled training data
#' # Using default KEGG file
#' training_results <- scCulture(
#'   data_dir = "./DATA_labeled",
#'   output_dir = "./training_results",
#'   mode = "build"
#' )
#'
#' # Or specify custom KEGG file
#' training_results_custom <- scCulture(
#'   data_dir = "./DATA_labeled",
#'   kegg_file = "custom_pathways.keg",
#'   output_dir = "./training_results",
#'   mode = "build"
#' )
#'
#' # Access training results
#' fingerprint_file <- training_results$fingerprint_file
#' training_accuracy <- training_results$evaluation_results$overall_accuracy
#' }
#'
#' \dontrun{
#' # Predict mode - apply fingerprints to new unlabeled data
#' prediction_results <- scCulture(
#'   data_dir = "./DATA_unlabeled",
#'   output_dir = "./prediction_results",
#'   mode = "predict",
#'   fingerprint_file = fingerprint_file
#' )
#'
#' # Access predictions
#' predictions <- prediction_results$seurat_object$classification_pred
#' confidence_scores <- prediction_results$seurat_object$prediction_confidence
#' }
#'
#' \dontrun{
#' # Advanced usage with progress tracking and parallel processing
#' # Using default KEGG file
#' results <- scCulture(
#'   data_dir = "./DATA_edit",
#'   output_dir = "./results",
#'   mode = "build",
#'   experiment_id = "my_experiment",
#'   progress = TRUE,
#'   parallel = TRUE,
#'   verbose = TRUE
#' )
#' }
#' @seealso
#' \code{\link{load_data}} for data loading
#' \code{\link{preprocess_data}} for data preprocessing
#' \code{\link{reduce_dimensions}} for dimensionality reduction
#' \code{\link{parse_kegg_keg}} for KEGG pathway analysis
#' \code{\link{predict_by_similarity}} for similarity-based prediction
#' \code{\link{predict_by_svm}} for SVM-based prediction
#' \code{\link{evaluate_predictions}} for result evaluation
#' \code{\link{create_evaluation_plots}} for visualization
#'
#' @export
scCulture <- function(data_dir, kegg_file = NULL, output_dir,
                      mode = "build", fingerprint_file = NULL,
                      experiment_id = "experiment",
                      perform_tsne = TRUE,
                      progress = FALSE, parallel = FALSE, n_cores = NULL,
                      verbose = TRUE) {
  # Input validation
  if (!is.character(data_dir) || length(data_dir) != 1) {
    stop("data_dir must be a single character string")
  }
  if (!dir.exists(data_dir)) {
    stop(sprintf("Data directory not found: %s", data_dir))
  }

  # Validate mode parameter
  if (!is.character(mode) || length(mode) != 1) {
    stop("mode must be a single character string")
  }
  if (!mode %in% c("build", "predict")) {
    stop("mode must be either 'build' or 'predict'")
  }

  # Mode-specific validation
  if (mode == "build") {
    # Use default KEGG file if none provided
    if (is.null(kegg_file)) {
      kegg_file <- system.file("extdata", "kegg", "sce00001.keg",
        package = "scCulturePredict"
      )
      if (!file.exists(kegg_file) || kegg_file == "") {
        stop("Default KEGG file not found. Please provide a kegg_file path.")
      }
      if (verbose) {
        message(sprintf("Using default KEGG file: %s", kegg_file))
      }
    }

    if (!is.character(kegg_file) || length(kegg_file) != 1) {
      stop("kegg_file must be a single character string for build mode")
    }
    if (!file.exists(kegg_file)) {
      stop(sprintf("KEGG file not found: %s", kegg_file))
    }
  } else if (mode == "predict") {
    if (is.null(fingerprint_file) || !is.character(fingerprint_file) || length(fingerprint_file) != 1) {
      stop("fingerprint_file must be provided as a single character string for predict mode")
    }
    if (!file.exists(fingerprint_file)) {
      stop(sprintf("Fingerprint file not found: %s", fingerprint_file))
    }
  }

  if (!is.character(output_dir) || length(output_dir) != 1) {
    stop("output_dir must be a single character string")
  }

  if (!is.character(experiment_id) || length(experiment_id) != 1) {
    stop("experiment_id must be a single character string")
  }

  if (!is.logical(perform_tsne) || length(perform_tsne) != 1) {
    stop("perform_tsne must be a single logical value")
  }

  if (!is.logical(progress) || length(progress) != 1) {
    stop("progress must be a single logical value")
  }

  if (!is.logical(parallel) || length(parallel) != 1) {
    stop("parallel must be a single logical value")
  }

  if (!is.null(n_cores)) {
    if (!is.numeric(n_cores) || length(n_cores) != 1 || n_cores < 1) {
      stop("n_cores must be a single positive number or NULL")
    }
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Set up parallel processing if requested
  if (parallel) {
    if (is.null(n_cores)) {
      n_cores <- parallel::detectCores() - 1
    }
    if (verbose) message(sprintf("Setting up parallel processing with %d cores...", n_cores))
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl))
  }

  # Set up progress bar if requested
  pb <- NULL
  if (progress) {
    pb <- txtProgressBar(min = 0, max = 100, style = 3, title = "scCulturePredict Analysis Progress")
    message("")
  }

  # Run pipeline with error handling
  tryCatch(
    {
      # Create output directory
      if (verbose) message("Creating output directory...")
      create_dir_if_not_exists(output_dir)
      if (progress) {
        setTxtProgressBar(pb, 5)
        message("")
      }

      # Load data (10%)
      if (verbose) message("Loading single-cell data...")
      seurat_object <- load_data(data_dir, experiment_id, verbose = verbose)
      if (progress) {
        setTxtProgressBar(pb, 15)
        message("")
      }

      # Preprocess data (25%)
      if (verbose) message("Preprocessing data...")
      seurat_object <- preprocess_data(seurat_object, verbose = verbose)
      if (progress) {
        setTxtProgressBar(pb, 25)
        message("")
      }

      # Perform dimensionality reduction (45%)
      if (verbose) message("Performing dimensionality reduction...")
      seurat_object <- reduce_dimensions(seurat_object, perform_tsne = perform_tsne, verbose = verbose)
      if (verbose) message("Dimensionality reduction complete")
      if (progress) {
        setTxtProgressBar(pb, 45)
        message("")
      }

      # Mode-specific pipeline execution
      if (mode == "build") {
        if (verbose) message("Running BUILD mode: generating fingerprints from labeled data...")
        result <- run_build_mode(seurat_object, kegg_file, output_dir, progress, pb, verbose)
      } else {
        if (verbose) message("Running PREDICT mode: applying fingerprints to unlabeled data...")
        result <- run_predict_mode(seurat_object, fingerprint_file, output_dir, progress, pb, verbose)
      }

      # Complete progress bar and return results
      if (progress) {
        if (verbose) message("Analysis complete!")
        setTxtProgressBar(pb, 100)
        close(pb)
        message("")
      }

      if (verbose) message("scCulturePredict analysis completed successfully!")

      return(result)
    },
    error = function(e) {
      if (progress && !is.null(pb)) {
        close(pb)
      }
      stop(sprintf("scCulturePredict analysis failed: %s", e$message))
    }
  )
}

#' Helper function for build mode
#' @keywords internal
#'
#' @return A list containing the complete results of BUILD mode analysis:
#'   \describe{
#'     \item{seurat_object}{Seurat object with predictions and evaluation metadata}
#'     \item{models}{List containing trained direct and SVM classification models}
#'     \item{evaluation_results}{List with accuracy metrics, confusion matrices, and performance statistics}
#'     \item{plots}{Named list of generated visualization plots}
#'     \item{fingerprints}{Matrix of calculated transcriptomic fingerprints for each culture condition}
#'     \item{pathway_activities}{Matrix of pathway activity scores used for training}
#'   }
#' @examples
#' \donttest{
#' # Example requires prepared data files
#' # result <- run_build_mode(data_dir = "path/to/data",
#' #                         experiment_id = "exp1",
#' #                         output_dir = tempdir())
#' }
run_build_mode <- function(seurat_object, kegg_file, output_dir, progress, pb, verbose) {
  # Parse KEGG pathways (55%)
  if (verbose) message("Parsing KEGG pathways...")
  kegg_pathways <- parse_kegg_keg(kegg_file, verbose = verbose)
  if (progress) {
    setTxtProgressBar(pb, 55)
    message("")
  }

  # Build fingerprints (65%)
  if (verbose) message("Building transcriptomic fingerprints...")
  pathway_results <- build_fingerprints(seurat_object, kegg_pathways, verbose = verbose)
  if (progress) {
    setTxtProgressBar(pb, 65)
    message("")
  }

  # Make predictions for model training (75%)
  if (verbose) message("Training prediction models...")
  similarity_results <- predict_by_similarity(
    pathway_results$pathway_matrix,
    pathway_results$signature_matrix,
    verbose = verbose
  )
  svm_results <- predict_by_svm(
    pathway_results$pathway_matrix,
    seurat_object,
    verbose = verbose
  )
  if (progress) {
    setTxtProgressBar(pb, 75)
    message("")
  }

  # Add predictions to Seurat object
  seurat_object$predicted_sample_1 <- similarity_results$predicted_direct
  seurat_object$predicted_sample_2 <- similarity_results$predicted_threshold

  # Use SVM predictions if available, otherwise fall back to similarity
  if (all(is.na(svm_results$predictions))) {
    seurat_object$classification_pred <- similarity_results$predicted_direct
    if (verbose) message("Using similarity-based predictions (SVM unavailable)")
  } else {
    seurat_object$classification_pred <- svm_results$predictions
  }

  # Evaluate model performance (85%)
  if (verbose) message("Evaluating model performance...")
  evaluation_results <- evaluate_predictions(seurat_object)
  if (progress) {
    setTxtProgressBar(pb, 85)
    message("")
  }

  # Save fingerprints and models for future predictions (90%)
  if (verbose) message("Saving fingerprints and models...")
  fingerprint_file <- file.path(output_dir, "scCulturePredict_fingerprints.rds")
  fingerprint_data <- list(
    kegg_pathways = kegg_pathways,
    pathway_results = pathway_results,
    svm_model = svm_results$svm_model,
    similarity_model = list(
      pathway_matrix = pathway_results$pathway_matrix,
      signature_matrix = pathway_results$signature_matrix
    ),
    metadata = list(
      created_date = Sys.Date(),
      n_pathways = length(kegg_pathways),
      n_samples = ncol(seurat_object),
      training_accuracy = evaluation_results$overall_accuracy
    )
  )
  save_object(fingerprint_data, fingerprint_file)
  if (progress) {
    setTxtProgressBar(pb, 90)
    message("")
  }

  # Create visualizations and save results (95%)
  if (verbose) message("Creating visualizations...")
  create_evaluation_plots(seurat_object, results_dir = output_dir, verbose = verbose)

  if (verbose) message("Saving build results...")
  save_object(seurat_object, file.path(output_dir, "seurat_object.rds"))
  save_object(pathway_results, file.path(output_dir, "pathway_results.rds"))
  save_object(
    list(similarity = similarity_results, svm = svm_results),
    file.path(output_dir, "prediction_results.rds")
  )
  save_object(evaluation_results, file.path(output_dir, "evaluation_results.rds"))
  if (progress) {
    setTxtProgressBar(pb, 95)
    message("")
  }

  return(list(
    seurat_object = seurat_object,
    pathway_results = pathway_results,
    prediction_results = list(similarity = similarity_results, svm = svm_results),
    #' @examples
    #' \donttest{
    #' # Example requires reference data and query data
    #' # result <- run_predict_mode(ref_data_dir = "path/to/ref",
    #' #                           query_data_dir = "path/to/query",
    #' #                           output_dir = tempdir())
    #' }
    evaluation_results = evaluation_results,
    fingerprint_file = fingerprint_file
  ))
}

#'
#' @return A list containing the complete results of PREDICT mode analysis:
#'   \describe{
#'     \item{seurat_object}{Query Seurat object with added prediction results and metadata}
#'     \item{predictions}{Character vector of predicted culture conditions for each cell}
#'     \item{confidence_scores}{Numeric vector of prediction confidence scores (0-1)}
#'     \item{plots}{Named list of generated visualization plots}
#'     \item{pathway_activities}{Matrix of calculated pathway activities used for predictions}
#'     \item{similarity_scores}{Matrix of similarity scores between cells and reference fingerprints}
#'   }
#' Helper function for predict mode
#' @keywords internal
run_predict_mode <- function(seurat_object, fingerprint_file, output_dir, progress, pb, verbose) {
  # Load fingerprints and models (55%)
  if (verbose) message("Loading fingerprints and trained models...")
  fingerprint_data <- readRDS(fingerprint_file)
  if (progress) {
    setTxtProgressBar(pb, 55)
    message("")
  }

  # Validate fingerprint data structure
  required_components <- c("kegg_pathways", "pathway_results", "svm_model", "similarity_model")
  missing_components <- setdiff(required_components, names(fingerprint_data))
  if (length(missing_components) > 0) {
    stop(sprintf(
      "Invalid fingerprint file. Missing components: %s",
      paste(missing_components, collapse = ", ")
    ))
  }

  # Calculate pathway activities for new data (65%)
  if (verbose) message("Calculating pathway activities for new data...")
  pathway_activities <- calculate_pathway_activities(
    seurat_object,
    fingerprint_data$kegg_pathways,
    verbose = verbose
  )

  # Use pre-built signature matrix from fingerprints
  pathway_results <- list(
    pathway_matrix = pathway_activities,
    signature_matrix = fingerprint_data$similarity_model$signature_matrix
  )
  if (progress) {
    setTxtProgressBar(pb, 65)
    message("")
  }

  # Make predictions using trained models (75%)
  if (verbose) message("Making predictions with trained models...")

  # Similarity-based predictions
  similarity_results <- predict_by_similarity(
    pathway_results$pathway_matrix,
    fingerprint_data$similarity_model$signature_matrix,
    verbose = verbose
  )

  # SVM predictions using pre-trained model
  # Prepare data frame with sanitized column names for SVM prediction
  prediction_data <- as.data.frame(t(pathway_results$pathway_matrix))
  colnames(prediction_data) <- make.names(colnames(prediction_data))

  # Feature matching safeguard: ensure prediction data has all features expected by the model
  svm_model <- fingerprint_data$svm_model
  svm_predictions <- tryCatch(
    {
      # Check if we can extract model terms (features)
      if (inherits(svm_model, "svm") && !is.null(svm_model$terms)) {
        # Get feature names from SVM model
        model_features <- attr(svm_model$terms, "term.labels")

        # Check for missing features
        missing_features <- setdiff(model_features, colnames(prediction_data))
        if (length(missing_features) > 0) {
          if (verbose) {
            message(sprintf("Note: Adding %d missing features with default values for SVM compatibility", length(missing_features)))
          }
          # Add missing columns with NA values
          for (feature in missing_features) {
            prediction_data[[feature]] <- NA
          }
        }
      }

      # Make predictions with complete feature set
      predict(svm_model, newdata = prediction_data)
    },
    error = function(e) {
      if (verbose) {
        message(sprintf("SVM prediction encountered an issue (%s). Proceeding with similarity-based predictions.", e$message))
      }
      return(rep(NA, nrow(prediction_data)))
    }
  )

  svm_results <- list(
    predictions = svm_predictions,
    model = fingerprint_data$svm_model
  )
  if (progress) {
    setTxtProgressBar(pb, 75)
    message("")
  }

  # Add predictions to Seurat object
  seurat_object$predicted_sample_1 <- similarity_results$predicted_direct
  seurat_object$predicted_sample_2 <- similarity_results$predicted_threshold

  # Use SVM predictions if available, otherwise fall back to similarity
  if (all(is.na(svm_results$predictions))) {
    seurat_object$classification_pred <- similarity_results$predicted_direct
    if (verbose) message("Using similarity-based predictions (SVM unavailable)")
  } else {
    seurat_object$classification_pred <- svm_results$predictions
  }

  # Calculate prediction confidence scores (85%)
  if (verbose) message("Calculating prediction confidence...")
  confidence_scores <- calculate_prediction_confidence(
    pathway_results$pathway_matrix,
    fingerprint_data$similarity_model$pathway_matrix
  )
  seurat_object$prediction_confidence <- confidence_scores

  evaluation_results <- list(
    prediction_confidence = confidence_scores,
    fingerprint_source = fingerprint_file,
    prediction_date = Sys.Date(),
    n_predictions = ncol(seurat_object)
  )
  if (progress) {
    setTxtProgressBar(pb, 85)
    message("")
  }

  # Create prediction visualizations (90%)
  if (verbose) message("Creating prediction visualizations...")
  create_prediction_plots(seurat_object, results_dir = output_dir, verbose = verbose)
  if (progress) {
    setTxtProgressBar(pb, 90)
    message("")
  }

  # Save prediction results (95%)
  if (verbose) message("Saving prediction results...")
  save_object(seurat_object, file.path(output_dir, "predicted_seurat_object.rds"))
  save_object(pathway_results, file.path(output_dir, "applied_pathway_results.rds"))
  save_object(
    list(similarity = similarity_results, svm = svm_results),
    file.path(output_dir, "prediction_results.rds")
  )
  save_object(evaluation_results, file.path(output_dir, "prediction_evaluation.rds"))
  if (progress) {
    setTxtProgressBar(pb, 95)
    message("")
  }

  return(list(
    seurat_object = seurat_object,
    pathway_results = pathway_results,
    prediction_results = list(similarity = similarity_results, svm = svm_results),
    evaluation_results = evaluation_results,
    fingerprint_source = fingerprint_file
  ))
}

#' Calculate prediction confidence scores
#' @keywords internal
#'
#' @return A numeric vector of confidence scores between 0 and 1 for each prediction. Higher values indicate more confident predictions.
calculate_prediction_confidence <- function(new_pathway_matrix, reference_pathway_matrix) {
  # Calculate correlation-based confidence scores
  correlations <- cor(new_pathway_matrix, reference_pathway_matrix, method = "pearson")

  # Safe max calculation that avoids warnings for empty/NA rows
  confidence_scores <- apply(correlations, 1, function(x) {
    non_na_values <- x[!is.na(x)]
    if (length(non_na_values) == 0) {
      return(0) # Default confidence if no valid correlations
    } else {
      return(max(non_na_values))
    }
  })

  # Convert to 0-1 scale and handle any remaining edge cases
  confidence_scores <- pmax(0, confidence_scores)
  confidence_scores[is.na(confidence_scores) | is.infinite(confidence_scores)] <- 0

  return(confidence_scores)
}

#'
#' @return A named list of ggplot2 objects with prediction visualizations tailored to the analysis type and available metadata.
#' Create prediction-specific visualizations
#' @keywords internal
create_prediction_plots <- function(seurat_object, results_dir, verbose = TRUE) {
  # This function would create prediction-specific plots
  # For now, we'll use the existing evaluation plots function
  # In a full implementation, this would create specialized plots for predictions
  create_evaluation_plots(seurat_object, results_dir = results_dir, verbose = verbose)
}
