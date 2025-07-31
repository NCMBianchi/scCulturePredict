# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.99.5] - 2025-07-31

### Fixed
- Fixed GitHub Actions workflow by removing base R packages (parallel, methods, stats, utils, tools) from BiocManager installation commands
- Split dependency installation into smaller, more manageable steps for better error tracking
- Added explicit package installation check before running BiocCheck
- Added `ask = FALSE` parameter to BiocManager::install() calls to prevent interactive prompts
- Improved workflow reliability with step-by-step installation and error handling

## [0.99.4] - 2025-07-31

### Fixed
- Fixed GitHub Actions workflow to properly install all package dependencies before running checks
- Improved dependency installation order to ensure BiocManager packages are available
- Enhanced workflow reliability for automated testing and validation

## [0.99.3] - 2025-07-31

### Changed
- Lowered R version requirement from 4.4.0 to 4.3.0 for broader compatibility with GitHub Actions and CI/CD environments

## [0.99.2] - 2025-07-31

### Fixed
- Fixed syntax error in `train_cell_type_classifier` function where `seq_len(min)(...)` was incorrectly parenthesized
- Fixed `vapply()` calls in `build_fingerprints` and `calculate_pathway_activities` functions by changing to `lapply()` for variable-length outputs
- Restored accidentally removed F1 score calculation block in `evaluate_cell_type_predictions` function

### Changed
- Improved code indentation compliance using styler package (reduced from 12% to 5% non-compliant lines)
- Achieved 0 ERRORS status in BiocCheck validation for GitHub Actions compatibility

## [0.99.1] - 2025-07-04

### Added
- BiocCheck output folder pattern to `.Rbuildignore` to prevent build errors

### Fixed
- Fixed all `1:n` patterns in examples, replaced with `seq_len()` for BiocCheck compliance
- Fixed code indentation issues using styler package (reduced from 11% to 7% non-compliant lines)
- Corrected SVM column name handling in `predict_by_svm()` function

### Changed
- Improved BiocCheck compliance, achieving 0 ERRORS status
- Updated documentation with properly formatted examples

### Removed
- Removed temporary helper scripts (`fix_coding_practices.R`, `fix_sapply_todos.R`)

## [0.99.0] - 2025-06-20

### Added

#### Core Functionality
- **Dual-mode main function `scCulture()`** with BUILD and PREDICT modes
- **BUILD Mode**: Generate transferable transcriptomic fingerprints from labeled training data
- **PREDICT Mode**: Apply pre-built fingerprints to unlabeled datasets for culture condition prediction
- **KEGG pathway analysis** for biological interpretation and feature extraction
- **Similarity-based prediction algorithm** using correlation-based matching
- **Support Vector Machine (SVM) prediction** with automatic feature handling
- **Prediction confidence scoring** [0;1] with edge case handling
- **Transferable fingerprint files** containing the trained model for culture media prediction, and pathway data

#### Data Integration
- Support for **CSV and 10X Genomics data formats** with and without appropriate headers
- **Seurat workflow integration**
- **Robust preprocessing pipeline** with normalization and scaling
- **Dimensionality reduction** using PCA, UMAP, and t-SNE
- **Metadata handling** with automatic sample detection

#### User Interface
- **Progress tracking** with minimal progress bars and detailed status messages
- **Parallel processing support** for computationally intensive tasks
- **Comprehensive input validation** with informative error messages
- **Verbose output options** for detailed analysis tracking
- **Visualization compatibility** with custom plotting function `plot_scCulture()`

#### Technical Infrastructure
- **Intelligent SVM fallback mechanism** - automatically switches to similarity-based predictions when SVM fails during prediction
- **Seurat v5 layer compatibility** with proper layer detection and handling
- **Memory-efficient processing** for large single-cell datasets
- **Cross-platform compatibility** (Windows, macOS, Linux)

### Fixed

- None in this initial release

### Changed

- None in this initial release

### Security

#### Input Validation
- **Parameter type checking** for all user inputs
- **File existence validation** before processing
- **Data format verification** to prevent processing errors
- **Memory usage monitoring** for large dataset handling

### Known Issues

#### Documentation Date Discrepancy
- **DESCRIPTION file date error**: Version 0.99.0 was erroneously dated as 2024-04-28 in the DESCRIPTION file instead of the correct date 2025-06-20. This has been corrected in version 0.99.1.

#### SVM Prediction Limitations
- **Feature mismatch scenarios**: SVM prediction may fail when there are significant differences between training and prediction datasets (e.g., different gene coverage, pathway availability)
- **Impact on functionality**: All functionality is maintained even when SVM encounters issues, ensuring reliable predictions

#### Technical Limitations
- **Cross-dataset compatibility**: Best results achieved when training and prediction datasets have similar gene coverage
- **KEGG pathway dependency**: Requires KEGG pathway definitions for optimal performance
- **Seurat version compatibility**: Optimized for Seurat v4+ with v5 layer compatibility

### Deprecated

- None in this initial release

### Removed

- None in this initial release
