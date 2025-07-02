#!/usr/bin/env Rscript

# Pre-commit hook for scUMAP package
# This script runs linting and styling checks on R files that are about to be committed
# To use this script as a pre-commit hook:
# 1. Copy or symlink this file to .git/hooks/pre-commit
# 2. Make it executable: chmod +x .git/hooks/pre-commit

# Get staged files
staged_files <- system("git diff --cached --name-only --diff-filter=ACM", intern = TRUE)
r_files <- staged_files[grepl("\\.R$", staged_files)]

if (length(r_files) == 0) {
  cat("No R files to check.\n")
  quit(status = 0)
}

# Check if required packages are installed
required_packages <- c("lintr", "styler")
missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

if (length(missing_packages) > 0) {
  cat("Missing required packages:", paste(missing_packages, collapse = ", "), "\n")
  cat("Please install them using: install.packages(c(\"", paste(missing_packages, collapse = "\", \""), "\"))\n", sep = "")
  quit(status = 1)
}

# Run linting checks
cat("Running linting checks...\n")
lint_issues <- list()
for (file in r_files) {
  cat("  Checking", file, "\n")
  issues <- lintr::lint(file)
  if (length(issues) > 0) {
    lint_issues[[file]] <- issues
    for (issue in issues) {
      cat(sprintf("    %s:%d:%d: %s: %s\n", 
                 file, 
                 issue$line_number, 
                 issue$column_number, 
                 issue$type, 
                 issue$message))
    }
  }
}

# Apply styling if needed
cat("Applying code styling...\n")
styled_files <- c()
for (file in r_files) {
  cat("  Styling", file, "\n")
  style_result <- styler::style_file(file, style = styler::tidyverse_style(indent_by = 2))
  if (style_result$changed) {
    styled_files <- c(styled_files, file)
    # Re-stage the file
    system(paste("git add", file))
  }
}

if (length(styled_files) > 0) {
  cat("Styled and re-staged files:", paste(styled_files, collapse = ", "), "\n")
}

# Decide whether to abort the commit
if (length(lint_issues) > 0) {
  cat("\nLinting issues found. Please fix them before committing.\n")
  cat("(If you want to bypass this check, use git commit --no-verify)\n")
  quit(status = 1)
} else {
  cat("\nAll checks passed.\n")
  quit(status = 0)
}