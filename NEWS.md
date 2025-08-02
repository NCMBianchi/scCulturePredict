# scCulturePredict 0.99.16 (2025-08-01)

## Minimal Workflow for Maximum Reliability

### GitHub Actions Workflow
* Created minimal, focused workflow to diagnose persistent setup failures
* Single OS testing (Ubuntu only) to reduce complexity
* Direct package installation without r-lib/actions/setup-r-dependencies
* Removed all non-essential features (multi-OS, coverage, caching)
* Focus on getting basic R CMD check and BiocCheck working
* If successful, will gradually add features back

# scCulturePredict 0.99.15 (2025-08-01)

## Return to Proven r-lib/actions Workflow

### GitHub Actions Workflow
* Reverted to simpler r-lib/actions workflow that worked in versions 0.99.8-11
* Removed complex biocthis workflow that was causing early setup failures
* Maintained BiocCheck integration for Bioconductor compliance
* Simplified dependency installation to single-pass approach
* Removed Docker containers and complex caching mechanisms
* Tests should now run successfully as in earlier versions

# scCulturePredict 0.99.14 (2025-08-01)

## Simplified Workflow for Better Reliability

### GitHub Actions Workflow
* Removed Docker container approach for Linux to fix setup failures
* Simplified workflow to use standard R setup across all platforms
* Updated to Bioconductor 3.19 (current release) from 3.18
* Maintained R 4.3 compatibility as established in version 0.99.3
* Fixed infrastructure issues preventing tests from running
* Ensured consistent setup process for Ubuntu, macOS, and Windows

# scCulturePredict 0.99.13 (2025-08-01)

## Workflow Debugging and Fixes

### GitHub Actions Workflow
* Debugging multi-OS workflow failures from version 0.99.12
* Investigating container availability and R version compatibility
* Addressing platform-specific dependency installation issues
* Working to resolve BiocManager and repository configuration problems

# scCulturePredict 0.99.12 (2025-08-01)

## Bioconductor-Specific CI/CD Updates

### GitHub Actions Workflow
* Switched to Bioconductor-specific workflow based on biocthis package template
* Added multi-OS testing support (Linux with Bioconductor Docker, macOS, Windows)
* Integrated proper BiocCheck execution into the workflow
* Uses Bioconductor package repositories and versioning
* Temporarily disabled code coverage (pending Codecov token setup)
* Maintains compatibility with Bioconductor 3.18 and R 4.3

# scCulturePredict 0.99.11 (2025-08-01)

## GitHub Actions Improvements

### CI/CD Workflow
* Rewrote GitHub Actions workflow using r-lib/actions best practices
* Replaced manual package installation with r-lib/actions/setup-r-dependencies
* Used r-lib/actions/check-r-package for standardized R CMD check
* Simplified BiocCheck execution with direct Rscript calls
* Improved workflow reliability and maintainability

# scCulturePredict 0.99.10 (2025-08-01)

## GitHub Actions Improvements

### CI/CD Workflow
* Added _R_CHECK_FORCE_SUGGESTS_=false environment variable to fix R CMD check error
* Fixed "Package suggested but not available: 'devtools'" error during package checking
* Set explicit R_LIBS_USER path for consistent package installation

# scCulturePredict 0.99.9 (2025-08-01)

## GitHub Actions Improvements

### CI/CD Workflow
* Fixed empty log issues by reverting from Rscript to R -e for better command execution
* Changed R CMD check to use rcmdcheck::rcmdcheck() directly to avoid devtools dependency
* Updated test runner to use testthat::test_local() instead of devtools::test()
* Ensured devtools is installed as a suggested package for development environments
* Reformatted system dependencies installation for better readability

# scCulturePredict 0.99.8 (2025-08-01)

## GitHub Actions Improvements

### CI/CD Workflow
* Added pandoc to system dependencies for vignette building
* Changed R -e to Rscript -e for more reliable command execution
* Removed force = TRUE parameter from BiocManager::install() calls
* Fixed package installation workflow to ensure all dependencies are properly installed

# scCulturePredict 0.99.7 (2025-07-31)

## GitHub Actions Improvements

### CI/CD Workflow
* Consolidated package installation steps to ensure devtools is installed properly
* Made R CMD check step more robust by using rcmdcheck as fallback
* Fixed workflow execution issues causing empty installation logs

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
