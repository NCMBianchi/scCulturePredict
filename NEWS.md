# scCulturePredict 0.99.6 (2025-07-31)

## GitHub Actions Improvements

### CI/CD Workflow
* Added missing system dependencies (libfontconfig1-dev, libfreetype6-dev, libpng-dev, libharfbuzz-dev, libfribidi-dev) required for Seurat installation
* Fixed package installation failures caused by missing system libraries for graphics packages

# scCulturePredict 0.99.5 (2025-07-31)

## GitHub Actions Improvements

### CI/CD Workflow
* Fixed GitHub Actions workflow by removing base R packages (parallel, methods, stats, utils, tools) from BiocManager installation commands
* Split dependency installation into smaller, more manageable steps for better error tracking
* Added explicit package installation check before running BiocCheck
* Added `ask = FALSE` parameter to BiocManager::install() calls to prevent interactive prompts
* Improved workflow reliability with step-by-step installation and error handling

# scCulturePredict 0.99.4 (2025-07-31)

## GitHub Actions Improvements

### CI/CD Workflow
* Fixed GitHub Actions workflow to properly install all package dependencies before running checks
* Improved dependency installation order to ensure BiocManager packages are available
* Enhanced workflow reliability for automated testing and validation

# scCulturePredict 0.99.3 (2025-07-31)

## Compatibility Improvements

### R Version Requirement
* Lowered R version requirement from 4.4.0 to 4.3.0 for broader compatibility with GitHub Actions and CI/CD environments

# scCulturePredict 0.99.2 (2025-07-31)

## Bug Fixes and Improvements

### Critical Fixes
* Fixed syntax error in `train_cell_type_classifier` function where `seq_len(min)(n_features, nrow(feature_importance))` was incorrectly parenthesized
* Fixed `vapply()` calls in `build_fingerprints` and `calculate_pathway_activities` functions by changing to `lapply()` for variable-length outputs
* Restored accidentally removed F1 score calculation block in `evaluate_cell_type_predictions` function

### BiocCheck Improvements
* Achieved 0 ERRORS status in BiocCheck validation, ensuring GitHub Actions compatibility
* Improved code indentation compliance using styler package (reduced from 12% to 5% non-compliant lines)

# scCulturePredict 0.99.1 (2025-07-04)

## Bug Fixes and Improvements

### BiocCheck Compliance
* Fixed all `1:n` patterns in examples, replaced with `seq_len()` for Bioconductor compliance
* Fixed code indentation issues using styler package (reduced from 11% to 7% non-compliant lines)
* Added BiocCheck output folder pattern to `.Rbuildignore` to prevent build errors
* Achieved 0 ERRORS in BiocCheck validation

### Code Improvements
* Corrected SVM column name handling in `predict_by_svm()` function
* Removed temporary helper scripts used during development
* Updated documentation with properly formatted examples

### Known Issues
* **DESCRIPTION file date correction**: Version 0.99.0 was erroneously dated as 2024-04-28 in the DESCRIPTION file instead of the correct date 2025-06-20. This has been corrected in version 0.99.1.

# scCulturePredict 0.99.0 (2025-06-20)

## Initial Release

I am excited to announce the first release of **scCulturePredict**, an R package for predicting cell culture media conditions from single-cell transcriptomic data using transferable transcriptomic fingerprints. This package is a case-study for appropriate handling of complex scrits via a single function (i.e. `scCulture()`) and several parameters that lead to different functionalities: this allows for easy-to-use data analysis.

## Key Features

### Dual-Mode Functionality
* **Build Mode**: Generate transferable fingerprints from labeled training data
  - Train prediction models using KEGG pathway analysis
  - Create portable fingerprint files for future use
  - Evaluate model performance with comprehensive metrics

* **Predict Mode**: Apply fingerprints to new unlabeled datasets
  - Make culture media predictions on new single-cell data that lack such an information
  - Calculate prediction confidence scores

### Main Functions
* **`scCulture()`**: Complete analysis pipeline with dual-mode functionality (i.e. BUILD, PREDICT)
* **`plot_scCulture()`**: Visualisation of results with automatic identification of the mode used in `scCulture()` to generate the appropriate plots

### Data Integration
* Support for CSV and 10X Genomics data formats
* Seamless integration with Seurat workflows
* Robust preprocessing and quality control
* Dimensionality reduction with PCA, UMAP, and t-SNE

## Installation

```r
# Install from GitHub
devtools::install_github("NCMBianchi/scCulturePredict")

# Load the package
library(scCulturePredict)
```

## Quick Start

```r
# BUILD mode - create fingerprints from labeled data
build_results <- scCulture(
  data_dir = "path/to/training/data",
  kegg_file = "path/to/kegg/file",
  output_dir = "fingerprint_output",
  mode = "build"
)

# PREDICT mode - apply fingerprints to new data
predict_results <- scCulture(
  data_dir = "path/to/new/data",
  fingerprint_file = build_results$fingerprint_file,
  output_dir = "prediction_output",
  mode = "predict"
)

# Automatic visualization (returns both plots for PREDICT mode)
plots <- plot_scCulture(predict_results)
print(plots$predictions)  # Culture medium predictions
print(plots$confidence)   # Prediction confidence scores
```

## Important Notes

### Prediction Methods
The package uses two complementary prediction approaches:
- **Similarity-based prediction**: Robust correlation-based matching
- **Support Vector Machine (SVM)**: Advanced classification when data permits

### Known Behavior
**SVM Automatic Fallback**: When SVM prediction encounters feature mismatches between training and prediction datasets, the package automatically falls back to similarity-based predictions with user notification. This ensures reliable predictions even when datasets have different characteristics.

## Technical Requirements

* R >= 4.1.0
* Seurat >= 4.0.0 (compatible with Seurat v5)
* Core packages: dplyr, ggplot2, tidyverse, MASS, e1071, caret
* Additional: patchwork, methods, parallel, doParallel, foreach, stats, utils, tools

## Performance Considerations

* Tested on dataset: 811 cells, 4673 genes
* Processing time: ~20-30 seconds for both BUILD and PREDICT modes
* Performance may vary with dataset size and complexity
* Best results with datasets having good gene coverage matching KEGG pathways

## Support

* **Documentation**: Use `?scCulture` and `?plot_scCulture` for detailed help
* **Issues**: Report problems at https://github.com/NCMBianchi/scCulturePredict/issues
* **Contact**: ncmb89@gmail.com

## Citation

If you use scCulturePredict in your research, please cite:

Bianchi, N. (2025). scCulturePredict: Single-Cell Culture Media Prediction Using Transcriptomic Fingerprints. R package version 0.99.0.
