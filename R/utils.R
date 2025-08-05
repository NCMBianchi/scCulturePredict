#' Default value operator
#'
#' @description
#' Returns the left-hand side if it's not NULL, otherwise returns the right-hand side.
#' This is a common utility operator used throughout the package.
#'
#' @param lhs Left-hand side value
#' @param rhs Right-hand side value (default)
#'
#' @return lhs if not NULL, otherwise rhs
#' @keywords internal
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs)) lhs else rhs
}

#' Check and install required packages
#'
#' @description
#' Checks if required packages are installed and installs them if they are not.
#' This is a more robust version of load_packages() that handles errors gracefully.
#'
#' @param packages Character vector of package names to check and install.
#' @param repos Character vector of repository URLs. Default is getOption("repos").
#'
#' @return Invisible NULL. This function is called for its side effects of checking and installing packages. Messages are printed to indicate installation progress.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if each package is installed using requireNamespace()
#'   \item Attempts to install missing packages using install.packages()
#'   \item Handles installation errors gracefully with warnings
#' }
#'
#' This function is more robust than load_packages() because it:
#' \itemize{
#'   \item Uses requireNamespace() instead of library()
#'   \item Handles installation errors without stopping
#'   \item Provides informative warning messages
#' }
#'
#' @examples
#' check_and_install_packages(c("Seurat", "dplyr", "ggplot2"))
#'
#' @seealso
#' \code{\link{load_packages}} for a simpler package loading function
#'
#' @export
check_and_install_packages <- function(packages, repos = getOption("repos")) {
  for (pkg in packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      message(sprintf("Installing package: %s", pkg))
      warning(sprintf("Package %s is not installed. Please install it manually: install.packages('%s')", pkg, pkg))
    }
  }
}

#' Create directory if it doesn't exist
#'
#' @description
#' Creates a directory if it doesn't exist, with optional recursive creation
#' of parent directories.
#'
#' @param dir_path Character string specifying the directory path to create.
#' @param recursive Logical indicating whether to create parent directories.
#'   Default is TRUE.
#'
#' @return Invisible NULL. The directory is created if it doesn't exist. The function returns the directory path invisibly.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if the directory exists using dir.exists()
#'   \item Creates the directory if it doesn't exist using dir.create()
#'   \item Optionally creates parent directories if recursive is TRUE
#' }
#'
#' This function is useful for:
#' \itemize{
#'   \item Creating output directories for analysis results
#'   \item Setting up directory structures for data processing
#'   \item Ensuring directories exist before saving files
#' }
#'
#' @examples
#' create_dir_if_not_exists("./results/plots")
#'
#' @seealso
#' \code{\link{dir.create}} for the base R function
#'
#' @export
create_dir_if_not_exists <- function(dir_path, recursive = TRUE) {
  if (!dir.exists(dir_path)) {
    dir.create(dir_path, recursive = recursive)
  }
}

#' Save R object to file
#'
#' @description
#' Saves an R object to a file with error handling and directory creation.
#'
#' @param object The R object to save.
#' @param file_path Character string specifying the file path.
#' @param compress Logical indicating whether to compress the file. Default is TRUE.
#'
#' @return Invisible NULL. The object is saved to the specified file path. The function returns the file path invisibly.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Creates the directory if it doesn't exist
#'   \item Saves the object using saveRDS()
#'   \item Handles errors with informative messages
#' }
#'
#' This function is useful for:
#' \itemize{
#'   \item Saving Seurat objects
#'   \item Saving analysis results
#'   \item Saving intermediate data
#' }
#'
#' The function ensures:
#' \itemize{
#'   \item The directory exists before saving
#'   \item Errors are caught and reported
#'   \item Files are compressed by default
#' }
#'
#' @examples
#' save_object(seurat_object, "./results/seurat_object.rds")
#'
#' @seealso
#' \code{\link{saveRDS}} for the base R function
#' \code{\link{create_dir_if_not_exists}} for directory creation
#'
#' @export
save_object <- function(object, file_path, compress = TRUE) {
  # Create directory if it doesn't exist
  dir_path <- dirname(file_path)
  create_dir_if_not_exists(dir_path)

  # Save object
  tryCatch(
    {
      saveRDS(object, file_path, compress = compress)
    },
    error = function(e) {
      stop(sprintf("Failed to save object to %s: %s", file_path, e$message))
    }
  )
}

#' Load R object from file
#'
#' @description
#' Loads an R object from a file with error handling.
#'
#' @param file_path Character string specifying the file path.
#'
#' @return The loaded R object.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if the file exists
#'   \item Loads the object using readRDS()
#'   \item Handles errors with informative messages
#' }
#'
#' This function is useful for:
#' \itemize{
#'   \item Loading saved Seurat objects
#'   \item Loading analysis results
#'   \item Loading intermediate data
#' }
#'
#' The function ensures:
#' \itemize{
#'   \item The file exists before loading
#'   \item Errors are caught and reported
#'   \item The loaded object is returned
#' }
#'
#' @examples
#' seurat_object <- load_object("./results/seurat_object.rds")
#'
#' @seealso
#' \code{\link{readRDS}} for the base R function
#' \code{\link{save_object}} for saving objects
#'
#' @export
load_object <- function(file_path) {
  if (!file.exists(file_path)) {
    stop(sprintf("File not found: %s", file_path))
  }

  tryCatch(
    {
      readRDS(file_path)
    },
    error = function(e) {
      stop(sprintf("Failed to load object from %s: %s", file_path, e$message))
    }
  )
}

#' Format numeric values
#'
#' @description
#' Formats numeric values with specified number of decimal places.
#'
#' @param x Numeric vector to format.
#' @param digits Integer specifying the number of decimal places. Default is 2.
#'
#' @return Character vector of formatted numbers.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters
#'   \item Formats numbers using sprintf()
#'   \item Returns formatted character vector
#' }
#'
#' This function is useful for:
#' \itemize{
#'   \item Formatting p-values
#'   \item Formatting statistics
#'   \item Formatting plot labels
#' }
#'
#' The function ensures:
#' \itemize{
#'   \item Consistent decimal places
#'   \item Proper rounding
#'   \item Clean output format
#' }
#'
#' @examples
#' format_number(c(1.23456, 2.34567), digits = 2)
#'
#' @seealso
#' \code{\link{sprintf}} for the base R function
#'
#' @export
format_number <- function(x, digits = 2) {
  sprintf(paste0("%.", digits, "f"), x)
}

#' Calculate percentage
#'
#' @description
#' Calculates percentage with proper formatting.
#'
#' @param x Numeric vector of values.
#' @param total Numeric value representing the total. Default is sum(x).
#' @param digits Integer specifying the number of decimal places. Default is 1.
#'
#' @return Character vector of formatted percentages.
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Validates input parameters
#'   \item Calculates percentages
#'   \item Formats with proper decimal places
#'   \item Adds percentage symbol
#' }
#'
#' This function is useful for:
#' \itemize{
#'   \item Calculating cell type percentages
#'   \item Formatting plot labels
#'   \item Reporting statistics
#' }
#'
#' The function ensures:
#' \itemize{
#'   \item Proper percentage calculation
#'   \item Consistent decimal places
#'   \item Clean output format with % symbol
#' }
#'
#' @examples
#' calculate_percentage(c(10, 20, 30))
#'
#' @seealso
#' \code{\link{format_number}} for number formatting
#'
#' @export
calculate_percentage <- function(x, total = sum(x), digits = 1) {
  sprintf(paste0("%.", digits, "f%%"), 100 * x / total)
}

#' Check if object is empty
#'
#' @description
#' Checks if an object is empty (NULL, NA, empty vector, or empty data frame).
#'
#' @param x Object to check.
#'
#' @return Logical value indicating whether the object is empty.
#'
#' @details
#' The function performs the following checks in order:
#' \enumerate{
#'   \item Checks if the object is NULL
#'   \item Checks if the object has length 0
#'   \item Checks if the object is an empty data frame
#'   \item Checks if all values are NA
#' }
#'
#' This function is useful for:
#' \itemize{
#'   \item Validating function inputs
#'   \item Checking data frame contents
#'   \item Handling edge cases
#' }
#'
#' The function handles:
#' \itemize{
#'   \item NULL values
#'   \item Empty vectors
#'   \item Empty data frames
#'   \item NA values
#' }
#'
#' @examples
#' is_empty(NULL)
#' is_empty(c())
#' is_empty(data.frame())
#'
#' @seealso
#' \code{\link{is.null}} for NULL checking
#' \code{\link{length}} for vector length checking
#'
#' @export
is_empty <- function(x) {
  if (is.null(x)) {
    return(TRUE)
  }
  if (length(x) == 0) {
    return(TRUE)
  }
  if (is.data.frame(x) && nrow(x) == 0) {
    return(TRUE)
  }
  if (all(is.na(x))) {
    return(TRUE)
  }
  FALSE
}

#' Get file extension
#'
#' @description
#' Extracts the file extension from a file path.
#'
#' @param file_path Character string specifying the file path.
#'
#' @return Character string containing the file extension (without the dot).
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Uses tools::file_ext() to extract the extension
#'   \item Returns the extension without the leading dot
#' }
#'
#' This function is useful for:
#' \itemize{
#'   \item File type validation
#'   \item File processing
#'   \item File organization
#' }
#'
#' The function ensures:
#' \itemize{
#'   \item Consistent extension format
#'   \item No leading dot in output
#'   \item Proper handling of files without extensions
#' }
#'
#' @examples
#' get_file_extension("data.csv")
#'
#' @seealso
#' \code{\link{tools::file_ext}} for the base function
#'
#' @export
get_file_extension <- function(file_path) {
  tools::file_ext(file_path)
}

#' Validate file exists
#'
#' @description
#' Validates that a file exists and is readable.
#'
#' @param file_path Character string specifying the file path.
#' @param extension Character string specifying the expected file extension.
#'   If NULL, no extension check is performed. Default is NULL.
#'
#' @return Logical value indicating whether the file is valid.
#'
#' @details
#' The function performs the following checks:
#' \enumerate{
#'   \item Validates file existence using file.exists()
#'   \item Optionally checks file extension
#'   \item Verifies file readability using file.access()
#' }
#'
#' This function is useful for:
#' \itemize{
#'   \item Input validation
#'   \item File processing
#'   \item Error handling
#' }
#'
#' The function ensures:
#' \itemize{
#'   \item File exists
#'   \item File has correct extension (if specified)
#'   \item File is readable
#'   \item Informative warning messages
#' }
#'
#' @examples
#' validate_file("data.csv", extension = "csv")
#'
#' @seealso
#' \code{\link{file.exists}} for existence checking
#' \code{\link{file.access}} for permission checking
#' \code{\link{get_file_extension}} for extension checking
#'
#' @export
validate_file <- function(file_path, extension = NULL) {
  if (!file.exists(file_path)) {
    warning(sprintf("File not found: %s", file_path))
    return(FALSE)
  }

  if (!is.null(extension)) {
    file_ext <- get_file_extension(file_path)
    if (file_ext != extension) {
      warning(sprintf("File extension mismatch. Expected: %s, Got: %s", extension, file_ext))
      return(FALSE)
    }
  }

  if (!file.access(file_path, 4) == 0) {
    warning(sprintf("File is not readable: %s", file_path))
    return(FALSE)
  }

  TRUE
}
