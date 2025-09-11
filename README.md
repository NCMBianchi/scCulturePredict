# scCulturePredict

[![Bioconductor](https://img.shields.io/badge/Bioconductor-under%20review-yellow)](https://github.com/Bioconductor/Contributions/issues)
[![Version](https://img.shields.io/badge/Version-0.99.32-orange)](https://github.com/nccb/scCulturePredict)
[![R-CMD-check-BiocCheck](https://github.com/NCMBianchi/scCulturePredict/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/NCMBianchi/scCulturePredict/actions/workflows/check-bioc.yml)
[![codecov](https://codecov.io/gh/NCMBianchi/scCulturePredict/branch/main/graph/badge.svg)](https://codecov.io/gh/NCMBianchi/scCulturePredict)
[![R](https://img.shields.io/badge/R-%3E%3D4.1.0-blue)](https://www.r-project.org/)
[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Build and Apply Transcriptomic Fingerprints for Cell Classification

scCulturePredict is an R package that provides dual functionality for classifying cells based on metabolic pathway signatures from single-cell transcriptomic data. While originally designed for culture media prediction, **scCulturePredict can classify cells based on any discrete metadata variable** (e.g., cell type, treatment condition, disease state, donor, timepoint, etc.) using metabolic pathway signatures.

**BUILD mode** generates transferable transcriptomic fingerprints from labeled training data, while **PREDICT mode** applies these pre-built fingerprints to unlabeled datasets for classification.

![](man/figures/scCulturePredict-banner.png)

## Features

### BUILD Mode (Generate Fingerprints)
- Train on labeled single-cell datasets
- Generate transferable transcriptomic fingerprints using KEGG pathway analysis
- Train both similarity-based and SVM prediction models
- Evaluate model performance with cross-validation
- Save fingerprints and models for future predictions

### PREDICT Mode (Apply Fingerprints)
- Apply pre-built fingerprints to unlabeled datasets
- Make culture media predictions using trained models
- Calculate prediction confidence scores
- Generate prediction-specific visualizations

### Core Capabilities
- Load and preprocess single-cell data (10X Genomics format or [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment) objects)
- Perform dimensionality reduction with UMAP and t-SNE
- Integrate with both [Seurat](https://satijalab.org/seurat/) and [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment) workflows
- Cross-dataset prediction with flexible pathway matching
- Comprehensive evaluation and visualization tools

![](man/figures/scCulturePredict-screenshot.png)

## Installation

### From Bioconductor (currently under review)

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("scCulturePredict")
```

### From GitHub (development version)

```r
# install.packages("devtools")
devtools::install_github("nccb/scCulturePredict")
```

## Quick Start

### BUILD Mode: Generate Fingerprints from Labeled Data

```r
library(scCulturePredict)

# Build fingerprints from labeled training data
# The function accepts both 10X Genomics data (via data_dir) and
# SingleCellExperiment objects (via sce_object)

# Option 1: Using 10X Genomics data
training_results <- scCulture(
  tenx_data_dir = "./DATA_labeled",      # Path to 10X data directory
  input_type = "10x",
  kegg_file = "kegg_pathways.keg",
  output_dir = "./training_results",
  mode = "build",
  experiment_id = "training",
  progress = TRUE,
  verbose = TRUE
)

# Option 2: Using SingleCellExperiment data
# Note: sce_object can be either a path to an RDS file OR an actual SCE object
# training_results <- scCulture(
#   sce_data_path = "labeled_cells.rds",  # Path to RDS file containing SCE object
#   input_type = "sce",
#   kegg_file = "kegg_pathways.keg",
#   output_dir = "./training_results",
#   mode = "build",
#   experiment_id = "training"
# )

# Access training results
fingerprint_file <- training_results$fingerprint_file
training_accuracy <- training_results$evaluation_results$overall_accuracy
print(paste("Training accuracy:", training_accuracy))
```

### PREDICT Mode: Apply Fingerprints to New Data

```r
# Apply fingerprints to unlabeled data
# Works with both 10X Genomics data and SingleCellExperiment objects

# Option 1: Using 10X Genomics data
prediction_results <- scCulture(
  tenx_data_dir = "./DATA_unlabeled",      # Path to 10X data directory
  input_type = "10x",
  output_dir = "./prediction_results",
  mode = "predict",
  fingerprint_file = fingerprint_file,  # From BUILD mode
  experiment_id = "predictions"
)

# Option 2: Using SingleCellExperiment data
# Note: sce_object can be either a path to an RDS file OR an actual SCE object
# prediction_results <- scCulture(
#   sce_data_path = "unlabeled_cells.rds",  # Path to RDS file containing SCE object
#   input_type = "sce",
#   output_dir = "./prediction_results",
#   mode = "predict",
#   fingerprint_file = fingerprint_file,
#   experiment_id = "predictions"
# )

# Access predictions
predictions <- prediction_results$seurat_object$classification_pred
confidence_scores <- prediction_results$seurat_object$prediction_confidence

# View results
head(data.frame(
  cell_barcode = colnames(prediction_results$seurat_object),
  predicted_class = predictions,
  confidence = confidence_scores
))
```

### Flexible Classification: Beyond Culture Media

scCulturePredict can classify cells based on **any discrete metadata variable**. Here are examples for different classification tasks:

```r
# Example 1: Cell Type Classification
# Your SCE object has a "cell_type" column with T cells, B cells, NK cells, etc.
cell_type_fingerprints <- scCulture(
  sce_data_path = "pbmc_data.rds",  # Can be path to RDS file or actual SCE object
  input_type = "sce",
  kegg_file = "human_kegg.keg",
  output_dir = "./cell_type_analysis",
  mode = "build",
  sample_column = "cell_type",  # Specify which metadata column to use
  experiment_id = "cell_type_classification"
)

# Example 2: Treatment Response Classification
# Your data has "treatment" column with Control, DrugA, DrugB
treatment_fingerprints <- scCulture(
  tenx_data_dir = "./treatment_data",
  input_type = "10x",
  kegg_file = "kegg_pathways.keg",
  output_dir = "./treatment_analysis",
  mode = "build",
  sample_column = "treatment",
  experiment_id = "drug_response"
)

# Example 3: Disease State Classification
# Your data has "condition" column with Healthy, Mild, Severe
disease_fingerprints <- scCulture(
  sce_data_path = "patient_data.rds",  # Can be path to RDS file or actual SCE object
  input_type = "sce",
  kegg_file = "kegg_pathways.keg",
  output_dir = "./disease_analysis",
  mode = "build",
  sample_column = "condition",
  experiment_id = "disease_state"
)
```

### Complete Workflow Example

```r
# Step 1: Build fingerprints (training phase)
training_results <- scCulture(
  tenx_data_dir = "./DATA_labeled",
  input_type = "10x",
  kegg_file = "sce00001.keg",
  output_dir = "./results/training",
  mode = "build"
)

# Step 2: Apply to new data (prediction phase)
prediction_results <- scCulture(
  tenx_data_dir = "./DATA_unlabeled",
  input_type = "10x",
  output_dir = "./results/predictions",
  mode = "predict",
  fingerprint_file = training_results$fingerprint_file
)

# Check prediction confidence
summary(prediction_results$seurat_object$prediction_confidence)
table(prediction_results$seurat_object$classification_pred)
```



## Data Format Requirements

scCulturePredict requires single-cell RNA-seq data in **10X Genomics format**:
- `matrix.mtx.gz` or `matrix.mtx` - Gene expression matrix
- `barcodes.tsv.gz` or `barcodes.tsv` - Cell barcodes
- `features.tsv.gz` or `features.tsv` - Gene information
- `metadata.tsv.gz` or `metadata.tsv` (optional) - Cell metadata with sample information

### Preprocessing GSE165686 Data (Optional)

If you're working with GSE165686 format files that have malformed headers (e.g., "x" in the first row), you can use the included shell script to preprocess the data:

```bash
# Location: inst/scripts/transform_files.sh
# Make script executable
chmod +x inst/scripts/transform_files.sh

# Run preprocessing
./inst/scripts/transform_files.sh input_directory output_directory
```

This script will:
- Remove malformed "x" headers from barcodes and features files
- Rename GSE165686-formatted files to standard 10X names
- Handle gzip compression automatically

### Handling Duplicate Gene Names in SingleCellExperiment Data

When working with SingleCellExperiment objects, duplicate gene names may occur due to various reasons (e.g., multiple transcripts, isoforms, or data processing artifacts). The `scCulture()` function provides flexible options for handling duplicates through the `handle_duplicates` parameter:

```r
# Example: Handle duplicate genes when using SingleCellExperiment data
results <- scCulture(
  sce_data_path = "data_with_duplicates.rds",
  input_type = "sce",
  kegg_file = "kegg_pathways.keg",
  output_dir = "./results",
  mode = "build",
  handle_duplicates = "make_unique"  # Default behavior
)
```

Available options for `handle_duplicates`:
- **"make_unique"** (default): Appends .1, .2, etc. to duplicate gene names
- **"aggregate"**: Sums expression values for duplicate genes
- **"first"**: Keeps only the first occurrence of duplicate genes
- **"error"**: Stops with an informative error if duplicates are found

The function will issue a warning when duplicates are detected and handled:
```
Warning: Found 5 duplicate gene names. Handling with method: make_unique
Example duplicates: GENE1, GENE2, GENE3...
```

This parameter ensures robust processing of real-world datasets while maintaining flexibility for different use cases.

## Documentation

Comprehensive documentation is available in the package:

- `vignette("scCulturePredict-introduction")` - Introduction to scCulturePredict
- `vignette("scCulturePredict-visualization")` - Visualisation guide



## Code Quality

scCulturePredict implements several code quality measures to ensure robustness and maintainability:

### Linting

The package uses `lintr` for static code analysis. To run linting checks:

```r
# Install lintr if needed
# install.packages("lintr")

# Run linting on the package
lintr::lint_package()
```

A `.lintr` configuration file is included in the package root.

### Code Formatting

Code formatting follows the Bioconductor style guidelines and is enforced using `styler`:

```r
# Install styler if needed
# install.packages("styler")

# Apply styling to the package
styler::style_pkg(style = styler::tidyverse_style(indent_by = 2))
```

### Comprehensive Checks

Run the comprehensive check script to ensure the package is ready for Bioconductor submission:

```r
# From the package root directory
Rscript scripts/check_package.R
```

This will run:
- R CMD check (with --as-cran flag)
- BiocCheck
- Linting checks
- Test coverage analysis
- Vignette building
- Example code execution

## Development

### Pre-commit Hook

To enforce code quality during development, you can install the pre-commit hook:

```bash
# From the package root directory
cp scripts/pre-commit-hook.R .git/hooks/pre-commit
chmod +x .git/hooks/pre-commit
```

### Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

1. Fork the repository
2. Create your feature branch (`git checkout -b feature/amazing-feature`)
3. Run the code quality checks (`Rscript scripts/check_package.R`)
4. Commit your changes (`git commit -m 'Add some amazing feature'`)
5. Push to the branch (`git push origin feature/amazing-feature`)
6. Open a Pull Request

## Citation

If you use scCulturePredict in your research, please cite (bibtex format):

```bibtex
@Manual{scCulturePredict2025,
  title = {scCulturePredict: Single-Cell Feature Prediction Using Transcriptomic Fingerprints},
  author = {Niccolò Bianchi},
  year = {2025},
  note = {R package version 0.99.32},
  url = {https://github.com/ncmbianchi/scCulturePredict},
}
```

## License

This project is licensed under the MIT License - see the LICENSE file for details.
