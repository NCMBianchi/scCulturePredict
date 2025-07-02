#!/usr/bin/env Rscript

# Comprehensive check script for scCulturePredict package
# This script runs multiple checks to ensure the package is ready for Bioconductor submission
# Usage: Rscript scripts/check_package.R

# Function to check if a package is installed
is_package_installed <- function(pkg) {
  return(requireNamespace(pkg, quietly = TRUE))
}

# Function to install missing packages
install_if_missing <- function(pkg, bioc = FALSE) {
  if (!is_package_installed(pkg)) {
    cat(sprintf("Installing package: %s\n", pkg))
    if (bioc) {
      if (!is_package_installed("BiocManager")) {
        install.packages("BiocManager")
      }
      BiocManager::install(pkg)
    } else {
      install.packages(pkg)
    }
  }
}

# Ensure required packages are installed
required_packages <- c(
  "devtools", "roxygen2", "testthat", "knitr",
  "rmarkdown", "covr", "lintr", "styler", "rcmdcheck"
)
bioc_packages <- c("BiocCheck", "BiocStyle")

for (pkg in required_packages) {
  install_if_missing(pkg)
}

for (pkg in bioc_packages) {
  install_if_missing(pkg, bioc = TRUE)
}

# Load necessary libraries
library(devtools)
library(roxygen2)
library(covr)
library(rcmdcheck)

# Set working directory to package root
pkg_path <- "."
pkg_name <- "scUMAP"

# Define colors for output
RED <- "\033[0;31m"
GREEN <- "\033[0;32m"
YELLOW <- "\033[0;33m"
BLUE <- "\033[0;34m"
RESET <- "\033[0m"

print_header <- function(text) {
  cat("\n", BLUE, "==================================================", RESET, "\n")
  cat(BLUE, text, RESET, "\n")
  cat(BLUE, "==================================================", RESET, "\n\n")
}

print_success <- function(text) {
  cat(GREEN, "✓ ", text, RESET, "\n")
}

print_warning <- function(text) {
  cat(YELLOW, "⚠ ", text, RESET, "\n")
}

print_error <- function(text) {
  cat(RED, "✗ ", text, RESET, "\n")
}

# Check if any processes failed
any_failures <- FALSE

# Step 1: Update documentation
print_header("Updating documentation")
tryCatch(
  {
    devtools::document(pkg = pkg_path)
    print_success("Documentation updated successfully")
  },
  error = function(e) {
    print_error(paste("Failed to update documentation:", e$message))
    any_failures <- TRUE
  }
)

# Step 2: Run R CMD check
print_header("Running R CMD check (with --as-cran flag)")
check_results <- tryCatch(
  {
    rcmdcheck::rcmdcheck(
      path = pkg_path,
      args = c("--as-cran", "--no-manual"),
      error_on = "never"
    )

    # Process results
    if (length(check_results$errors) > 0) {
      for (err in check_results$errors) {
        print_error(paste("Error:", err))
      }
      any_failures <- TRUE
    } else {
      print_success("No errors in R CMD check")
    }

    if (length(check_results$warnings) > 0) {
      for (warn in check_results$warnings) {
        print_warning(paste("Warning:", warn))
      }
    } else {
      print_success("No warnings in R CMD check")
    }

    if (length(check_results$notes) > 0) {
      for (note in check_results$notes) {
        cat(YELLOW, "Note:", note, RESET, "\n")
      }
    } else {
      print_success("No notes in R CMD check")
    }

    check_results
  },
  error = function(e) {
    print_error(paste("R CMD check failed:", e$message))
    any_failures <- TRUE
    NULL
  }
)

# Step 3: Run BiocCheck
print_header("Running BiocCheck")
bioc_check_results <- tryCatch(
  {
    BiocCheck::BiocCheck(
      package = pkg_path,
      new_package = TRUE,
      build_output = TRUE,
      check_version = TRUE,
      force = TRUE
    )

    # Process results
    if (length(bioc_check_results$error) > 0) {
      for (err in bioc_check_results$error) {
        print_error(paste("BiocCheck Error:", err))
      }
      any_failures <- TRUE
    } else {
      print_success("No errors in BiocCheck")
    }

    if (length(bioc_check_results$warning) > 0) {
      for (warn in bioc_check_results$warning) {
        print_warning(paste("BiocCheck Warning:", warn))
      }
    } else {
      print_success("No warnings in BiocCheck")
    }

    if (length(bioc_check_results$note) > 0) {
      for (note in bioc_check_results$note) {
        cat(YELLOW, "BiocCheck Note:", note, RESET, "\n")
      }
    } else {
      print_success("No notes in BiocCheck")
    }

    bioc_check_results
  },
  error = function(e) {
    print_error(paste("BiocCheck failed:", e$message))
    any_failures <- TRUE
    NULL
  }
)

# Step 4: Run linting
print_header("Running linting checks")
r_files <- list.files(
  path = file.path(pkg_path, "R"),
  pattern = "\\.R$",
  full.names = TRUE
)

lint_issues <- list()
for (file in r_files) {
  issues <- lintr::lint(file)
  if (length(issues) > 0) {
    lint_issues[[file]] <- issues
    print_warning(sprintf("Linting issues in %s:", basename(file)))
    for (issue in issues) {
      cat(YELLOW, sprintf(
        "  Line %d, Col %d: %s - %s\n",
        issue$line_number,
        issue$column_number,
        issue$type,
        issue$message
      ), RESET)
    }
  }
}

if (length(lint_issues) == 0) {
  print_success("No linting issues found")
}

# Step 5: Run tests and check coverage
print_header("Running tests and checking coverage")
test_results <- tryCatch(
  {
    devtools::test(pkg = pkg_path)
    print_success("All tests passed")

    # Get test coverage
    coverage <- covr::package_coverage(path = pkg_path)
    cov_percent <- covr::percent_coverage(coverage)

    if (cov_percent < 80) {
      print_warning(sprintf("Test coverage is %.2f%%, which is below 80%%", cov_percent))
    } else {
      print_success(sprintf("Test coverage is %.2f%%", cov_percent))
    }

    # Show functions with low coverage
    zero_cov <- covr::zero_coverage(coverage)
    if (nrow(zero_cov) > 0) {
      print_warning("Functions with zero coverage:")
      for (i in 1:min(nrow(zero_cov), 10)) {
        cat(YELLOW, sprintf("  %s:%d\n", zero_cov$filename[i], zero_cov$line[i]), RESET)
      }
      if (nrow(zero_cov) > 10) {
        cat(YELLOW, sprintf("  ... and %d more\n", nrow(zero_cov) - 10), RESET)
      }
    }

    TRUE
  },
  error = function(e) {
    print_error(paste("Test failures:", e$message))
    any_failures <- TRUE
    FALSE
  }
)

# Step 6: Check vignettes
print_header("Building vignettes")
vignette_results <- tryCatch(
  {
    devtools::build_vignettes(pkg = pkg_path)
    print_success("Vignettes built successfully")
    TRUE
  },
  error = function(e) {
    print_error(paste("Failed to build vignettes:", e$message))
    any_failures <- TRUE
    FALSE
  }
)

# Step 7: Check examples
print_header("Running examples")
example_results <- tryCatch(
  {
    devtools::run_examples(pkg = pkg_path)
    print_success("All examples run successfully")
    TRUE
  },
  error = function(e) {
    print_error(paste("Failed to run examples:", e$message))
    any_failures <- TRUE
    FALSE
  }
)

# Final summary
print_header("Check Summary")
if (any_failures) {
  print_error("Some checks failed. Please address the issues before submitting to Bioconductor.")
} else {
  if (length(check_results$warnings) > 0 ||
    length(bioc_check_results$warning) > 0 ||
    length(lint_issues) > 0) {
    print_warning("All critical checks passed, but there are warnings to address.")
  } else {
    print_success("All checks passed successfully! The package is ready for submission.")
  }
}

# Provide a reminder about submission requirements
cat("\n", BLUE, "Bioconductor Submission Checklist:", RESET, "\n")
cat("1. Ensure all version numbers are correct in DESCRIPTION\n")
cat("2. Update NEWS.md with changes in this version\n")
cat("3. Make sure all vignettes build correctly\n")
cat("4. Verify dependencies are correctly specified\n")
cat("5. Ensure biocViews in DESCRIPTION are appropriate\n")
cat("6. Submit package via https://github.com/Bioconductor/Contributions\n")

# Return status code based on failures
if (any_failures) {
  quit(status = 1)
} else {
  quit(status = 0)
}
