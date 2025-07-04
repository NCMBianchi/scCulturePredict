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
