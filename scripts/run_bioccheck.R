#!/usr/bin/env Rscript

# Script to run BiocCheck on the scUMAP package
# Usage: Rscript scripts/run_bioccheck.R

# Ensure BiocCheck is installed
if (!requireNamespace("BiocCheck", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("BiocCheck")
}

# Load the necessary packages
library(BiocCheck)

# Get the path to the package
pkg_path <- "."

# Run BiocCheck with recommended settings
cat("Running BiocCheck on the scCulturePredict package...\n")
check_results <- BiocCheck(
  package = pkg_path,
  new_package = TRUE,
  build_output = TRUE,
  check_version = TRUE,
  force = TRUE,
  source_control = "git"
)

# Process and display the results
cat("\n=======================================\n")
cat("BiocCheck Results Summary:\n")
cat("=======================================\n")

# Print the summary
if (length(check_results$error) > 0) {
  cat("ERRORS (must be fixed):\n")
  for (err in check_results$error) {
    cat("  - ", err, "\n")
  }
}

if (length(check_results$warning) > 0) {
  cat("\nWARNINGS (should be addressed):\n")
  for (warn in check_results$warning) {
    cat("  - ", warn, "\n")
  }
}

if (length(check_results$note) > 0) {
  cat("\nNOTES (consider addressing):\n")
  for (note in check_results$note) {
    cat("  - ", note, "\n")
  }
}

# Print overall status
cat("\n=======================================\n")
if (length(check_results$error) == 0) {
  if (length(check_results$warning) == 0) {
    cat("SUCCESS: No errors or warnings found!\n")
  } else {
    cat("PARTIAL SUCCESS: No errors found, but there are warnings to address.\n")
  }
} else {
  cat("FAILURE: Errors must be fixed before submitting to Bioconductor.\n")
}
cat("=======================================\n")

# Return status code based on presence of errors
if (length(check_results$error) > 0) {
  quit(status = 1)
} else {
  quit(status = 0)
}
