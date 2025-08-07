#' scCulturePredict: Single-Cell Culture Media Prediction Package
#'
#' @description
#' The scCulturePredict package provides dual functionality for analyzing single-cell RNA sequencing data
#' to predict cell culture media conditions. The main function \code{\link{scCulture}} operates
#' in two modes: "build" to generate transferable transcriptomic fingerprints from labeled
#' training data, and "predict" to apply these fingerprints to unlabeled datasets. Individual
#' functions are available for data loading, preprocessing, KEGG pathway analysis, and prediction
#' using both similarity-based and machine learning approaches.
#'
#' @importFrom methods hasMethod
#' @importFrom stats chisq.test cor p.adjust predict sd t.test var
#' @importFrom utils head install.packages installed.packages read.csv read.table setTxtProgressBar txtProgressBar write.table
#'
#' @details
#' The package offers two primary workflows:
#'
#' **Build Mode** (\code{mode = "build"}):
#' \itemize{
#'   \item Generate transferable fingerprints from labeled training data
#'   \item Train similarity-based and SVM prediction models
#'   \item Evaluate model performance with cross-validation
#'   \item Save fingerprints and models for future use
#' }
#'
#' **Predict Mode** (\code{mode = "predict"}):
#' \itemize{
#'   \item Apply pre-built fingerprints to unlabeled datasets
#'   \item Make culture media predictions using trained models
#'   \item Calculate prediction confidence scores
#'   \item Generate prediction-specific visualizations
#' }
#'
#' **Core Functionalities**:
#' \itemize{
#'   \item Data loading with support for CSV and 10X Genomics formats
#'   \item Preprocessing including normalization and dimensionality reduction
#'   \item KEGG pathway analysis and transcriptomic fingerprint generation
#'   \item Dual prediction approaches: similarity-based and SVM machine learning
#'   \item Comprehensive evaluation and visualization tools
#' }
#'
#' For most users, the recommended approach is to use the \code{\link{scCulture}} function
#' which handles both build and predict workflows with a single function call.
#'
#' @section Main Pipeline:
#' The \code{\link{scCulture}} function provides dual-mode functionality:
#'
#' **Build Mode Pipeline**:
#' \itemize{
#'   \item Loads and preprocesses labeled training data
#'   \item Performs dimensionality reduction (PCA, UMAP, optional t-SNE)
#'   \item Conducts KEGG pathway analysis and builds fingerprints
#'   \item Trains similarity-based and SVM prediction models
#'   \item Evaluates model performance with comprehensive metrics
#'   \item Saves transferable fingerprints and trained models
#'   \item Creates training evaluation visualizations
#' }
#'
#' **Predict Mode Pipeline**:
#' \itemize{
#'   \item Loads pre-built fingerprints and trained models
#'   \item Loads and preprocesses unlabeled target data
#'   \item Applies fingerprints to generate pathway profiles
#'   \item Makes predictions using trained models
#'   \item Calculates prediction confidence scores
#'   \item Creates prediction-specific visualizations
#'   \item Saves prediction results and confidence metrics
#' }
#'
#' @section Data Loading:
#' The \code{\link{load_data}} function provides robust data loading capabilities:
#' \itemize{
#'   \item Supports both shell script and R-based file preparation
#'   \item Automatic handling of malformed input files
#'   \item Seamless integration with Seurat objects
#' }
#'
#' @section Preprocessing:
#' The preprocessing functions (\code{\link{preprocess_data}}, \code{\link{reduce_dimensions}})
#' handle data normalization and dimensionality reduction.
#'
#' @section KEGG Analysis:
#' KEGG pathway analysis is performed using \code{\link{parse_kegg_keg}} and
#' \code{\link{build_fingerprints}} functions.
#'
#' @section Prediction:
#' Cell type/state prediction can be done using either similarity-based
#' (\code{\link{predict_by_similarity}}) or SVM-based (\code{\link{predict_by_svm}})
#' approaches.
#'
#' @section Evaluation:
#' Results can be evaluated and visualized using \code{\link{evaluate_predictions}}
#' and \code{\link{create_evaluation_plots}} functions.
#'
#' @author
#' Your Name <your.email@example.com>
#'
#' @seealso
#' \code{\link{scCulture}} for complete analysis pipeline (recommended)
#' \code{\link{plot_scCulture}} for visualization of results
#' \code{\link{load_data}} for data loading
#' \code{\link{preprocess_data}} for data preprocessing
#' \code{\link{reduce_dimensions}} for dimensionality reduction
#' \code{\link{parse_kegg_keg}} for KEGG pathway analysis
#' \code{\link{predict_by_similarity}} for similarity-based prediction
#' \code{\link{predict_by_svm}} for SVM-based prediction
#' \code{\link{evaluate_predictions}} for result evaluation
#' \code{\link{create_evaluation_plots}} for creating evaluation plots

#' @keywords internal
"_PACKAGE"
