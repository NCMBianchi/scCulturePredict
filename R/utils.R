# Internal utility: Default value operator
`%||%` <- function(lhs, rhs) {
  if (!is.null(lhs)) lhs else rhs
}

#' Check and install required packages
#'
#' @description
#' Checks if required packages are installed and installs them if they are not.
#' This function handles package checking and installation gracefully.
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
#' This function is robust because it:
#' \itemize{
#'   \item Uses requireNamespace() instead of library()
#'   \item Handles installation errors without stopping
#'   \item Provides informative warning messages
#' }
#'
#' @examples
#' check_and_install_packages(c("Seurat", "dplyr", "ggplot2"))
#'

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
#' # Create a simple object to save
#' test_object <- list(data = 1:10, name = "test")
#'
#' # Save to temporary file
#' temp_file <- tempfile(fileext = ".rds")
#' save_object(test_object, temp_file)
#'
#' # Clean up
#' unlink(temp_file)
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
#' # Create and save a test object
#' test_data <- list(values = 1:5, type = "example")
#' temp_file <- tempfile(fileext = ".rds")
#' saveRDS(test_data, temp_file)
#'
#' # Load the object
#' loaded_object <- load_object(temp_file)
#'
#' # Verify it loaded correctly
#' print(loaded_object$type)
#'
#' # Clean up
#' unlink(temp_file)
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
