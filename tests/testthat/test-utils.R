# Test utility functions

test_that("check_and_install_packages works correctly", {
  # Test with packages that are definitely installed (base R packages)
  expect_silent({
    check_and_install_packages(c("base", "stats", "utils"))
  })

  # Test with empty vector
  expect_silent({
    check_and_install_packages(character(0))
  })

  # Test with NULL input
  expect_silent({
    check_and_install_packages(NULL)
  })

  # Test that it returns invisible NULL
  result <- check_and_install_packages(c("base"))
  expect_null(result)

  # Test with a single package
  expect_silent({
    check_and_install_packages("methods")
  })
})

test_that("load_object handles RDS files correctly", {
  # Create a temporary RDS file with test data
  temp_file <- tempfile(fileext = ".rds")
  test_data <- list(
    a = 1:10,
    b = "test string",
    c = data.frame(x = 1:3, y = c("a", "b", "c"))
  )
  saveRDS(test_data, temp_file)

  # Test loading the file
  loaded_data <- load_object(temp_file)
  expect_identical(loaded_data, test_data)

  # Test loading again to ensure consistency
  loaded_data2 <- load_object(temp_file)
  expect_identical(loaded_data2, test_data)

  # Clean up
  unlink(temp_file)
})

test_that("load_object handles errors correctly", {
  # Test with non-existent file
  expect_error(
    load_object("nonexistent_file.rds"),
    regexp = "File not found|does not exist"
  )

  # Test with invalid file path
  expect_error(
    load_object("/invalid/path/to/file.rds"),
    regexp = "File not found|does not exist"
  )

  # Test with NULL input
  expect_error(
    load_object(NULL),
    regexp = "File not found|must be|invalid"
  )

  # Test with empty string
  expect_error(
    load_object(""),
    regexp = "File not found|does not exist|invalid"
  )
})

test_that("load_object handles different data types", {
  # Test with various R objects
  temp_file <- tempfile(fileext = ".rds")

  # Test with a matrix
  test_matrix <- matrix(1:12, nrow = 3)
  saveRDS(test_matrix, temp_file)
  loaded_matrix <- load_object(temp_file)
  expect_identical(loaded_matrix, test_matrix)

  # Test with a data frame
  test_df <- data.frame(
    col1 = 1:5,
    col2 = letters[1:5],
    col3 = c(TRUE, FALSE, TRUE, FALSE, TRUE)
  )
  saveRDS(test_df, temp_file)
  loaded_df <- load_object(temp_file)
  expect_identical(loaded_df, test_df)

  # Test with a single value
  test_value <- 42
  saveRDS(test_value, temp_file)
  loaded_value <- load_object(temp_file)
  expect_identical(loaded_value, test_value)

  # Test with a complex nested list
  test_complex <- list(
    numbers = 1:10,
    strings = c("a", "b", "c"),
    nested = list(
      inner = list(value = 100),
      data = data.frame(x = 1:3, y = 4:6)
    )
  )
  saveRDS(test_complex, temp_file)
  loaded_complex <- load_object(temp_file)
  expect_identical(loaded_complex, test_complex)

  # Clean up
  unlink(temp_file)
})

test_that("save_object works correctly", {
  # Create test data
  test_data <- list(x = 1:10, y = letters[1:10])

  # Test saving to a temporary file
  temp_file <- tempfile(fileext = ".rds")

  # Test saving
  expect_silent({
    save_object(test_data, temp_file)
  })

  # Verify the file was created
  expect_true(file.exists(temp_file))

  # Verify the data can be loaded back correctly
  loaded_data <- readRDS(temp_file)
  expect_identical(loaded_data, test_data)

  # Test overwriting existing file
  expect_silent({
    save_object(test_data, temp_file)
  })

  # Clean up
  unlink(temp_file)
})

test_that("save_object creates directory if needed", {
  # Create a path with a non-existent directory
  temp_dir <- tempfile("test_dir")
  temp_file <- file.path(temp_dir, "test.rds")
  test_data <- list(a = 1, b = 2)

  # save_object should create the directory
  expect_silent({
    save_object(test_data, temp_file)
  })

  # Verify directory and file were created
  expect_true(dir.exists(temp_dir))
  expect_true(file.exists(temp_file))

  # Verify data integrity
  loaded_data <- readRDS(temp_file)
  expect_identical(loaded_data, test_data)

  # Clean up
  unlink(temp_dir, recursive = TRUE)
})

test_that("create_dir_if_not_exists works correctly", {
  # Test creating a new directory
  temp_dir <- tempfile("test_create_dir")
  expect_false(dir.exists(temp_dir))

  # Create directory
  expect_silent({
    create_dir_if_not_exists(temp_dir)
  })
  expect_true(dir.exists(temp_dir))

  # Test with existing directory (should not error)
  expect_silent({
    create_dir_if_not_exists(temp_dir)
  })

  # Test creating another new directory
  temp_dir2 <- tempfile("test_create_dir2")
  expect_silent({
    create_dir_if_not_exists(temp_dir2)
  })
  expect_true(dir.exists(temp_dir2))

  # Clean up
  unlink(temp_dir, recursive = TRUE)
  unlink(temp_dir2, recursive = TRUE)
})

test_that("%||% operator works correctly", {
  # Test with NULL on left side
  expect_equal(NULL %||% "default", "default")
  expect_equal(NULL %||% 42, 42)
  expect_equal(NULL %||% list(a = 1), list(a = 1))

  # Test with non-NULL on left side
  expect_equal("value" %||% "default", "value")
  expect_equal(10 %||% 42, 10)
  expect_equal(list(a = 1) %||% list(b = 2), list(a = 1))

  # Test with both NULL
  expect_null(NULL %||% NULL)

  # Test with empty values (should return the empty value, not the default)
  expect_equal(character(0) %||% "default", character(0))
  expect_equal(numeric(0) %||% 42, numeric(0))
  expect_equal(list() %||% list(a = 1), list())
})

test_that("validate_and_fix_file handles normal files correctly", {
  # Create a temporary file with normal data
  temp_file <- tempfile(fileext = ".tsv")
  test_data <- data.frame(
    gene = c("Gene1", "Gene2", "Gene3"),
    value1 = c(1, 2, 3),
    value2 = c(4, 5, 6)
  )
  write.table(test_data, temp_file, sep = "\t", row.names = FALSE, quote = FALSE)

  # Test reading the file
  result <- scCulturePredict:::validate_and_fix_file(temp_file, sep = "\t", header = TRUE, verbose = FALSE)
  expect_equal(result, test_data)

  # Test with verbose = TRUE
  expect_message(
    scCulturePredict:::validate_and_fix_file(temp_file, sep = "\t", header = TRUE, verbose = TRUE),
    "Reading file"
  )

  # Clean up
  unlink(temp_file)
})

test_that("validate_and_fix_file handles malformed headers", {
  # Create a file with malformed header (starting with 'x')
  temp_file <- tempfile(fileext = ".tsv")
  writeLines(c(
    "x",
    "gene\tvalue1\tvalue2",
    "Gene1\t1\t4",
    "Gene2\t2\t5",
    "Gene3\t3\t6"
  ), temp_file)

  # Test that it skips the malformed header
  result <- scCulturePredict:::validate_and_fix_file(temp_file, sep = "\t", header = TRUE, verbose = FALSE)
  expected <- data.frame(
    gene = c("Gene1", "Gene2", "Gene3"),
    value1 = c(1, 2, 3),
    value2 = c(4, 5, 6)
  )
  expect_equal(result, expected)

  # Test with verbose to check message
  expect_message(
    scCulturePredict:::validate_and_fix_file(temp_file, sep = "\t", header = TRUE, verbose = TRUE),
    "malformed header"
  )

  # Clean up
  unlink(temp_file)
})

test_that("validate_and_fix_file handles input validation", {
  # Test with non-existent file
  expect_error(
    scCulturePredict:::validate_and_fix_file("nonexistent.txt", verbose = FALSE),
    "File not found"
  )

  # Test with invalid file_path type
  expect_error(
    scCulturePredict:::validate_and_fix_file(123, verbose = FALSE),
    "file_path must be a single character string"
  )

  # Test with NULL file_path
  expect_error(
    scCulturePredict:::validate_and_fix_file(NULL, verbose = FALSE),
    "file_path must be a single character string"
  )

  # Test with vector of file paths
  expect_error(
    scCulturePredict:::validate_and_fix_file(c("file1.txt", "file2.txt"), verbose = FALSE),
    "file_path must be a single character string"
  )

  # Create a valid temp file for other parameter tests
  temp_file <- tempfile(fileext = ".txt")
  writeLines("col1,col2\n1,2", temp_file)

  # Test invalid sep parameter
  expect_error(
    scCulturePredict:::validate_and_fix_file(temp_file, sep = 123, verbose = FALSE),
    "sep must be a single character string"
  )

  # Test invalid header parameter
  expect_error(
    scCulturePredict:::validate_and_fix_file(temp_file, sep = ",", header = "yes", verbose = FALSE),
    "header must be a single logical value"
  )

  # Test invalid verbose parameter
  expect_error(
    scCulturePredict:::validate_and_fix_file(temp_file, sep = ",", header = TRUE, verbose = "yes"),
    "verbose must be a single logical value"
  )

  # Clean up
  unlink(temp_file)
})

test_that("validate_and_fix_file handles different separators", {
  # Test with comma-separated file
  temp_csv <- tempfile(fileext = ".csv")
  test_data <- data.frame(
    col1 = c("A", "B", "C"),
    col2 = c(1, 2, 3)
  )
  write.csv(test_data, temp_csv, row.names = FALSE, quote = FALSE)

  result <- scCulturePredict:::validate_and_fix_file(temp_csv, sep = ",", header = TRUE, verbose = FALSE)
  expect_equal(result, test_data)

  # Test with tab-separated file
  temp_tsv <- tempfile(fileext = ".tsv")
  write.table(test_data, temp_tsv, sep = "\t", row.names = FALSE, quote = FALSE)

  result <- scCulturePredict:::validate_and_fix_file(temp_tsv, sep = "\t", header = TRUE, verbose = FALSE)
  expect_equal(result, test_data)

  # Clean up
  unlink(temp_csv)
  unlink(temp_tsv)
})

test_that("process_metadata adds metadata to Seurat object", {
  # Skip if Seurat is not available
  skip_if_not_installed("Seurat")
  library(Seurat)

  # Create a minimal Seurat object
  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  colnames(counts) <- paste0("Cell", 1:10)
  rownames(counts) <- paste0("Gene", 1:10)
  seurat_obj <- CreateSeuratObject(counts = counts)

  # Create metadata data frame
  metadata <- data.frame(
    barcode = paste0("Cell", 1:10),
    condition = rep(c("A", "B"), each = 5),
    batch = rep(c(1, 2), 5),
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$barcode

  # Create barcodes data frame
  barcodes <- data.frame(
    barcode = paste0("Cell", 1:10),
    stringsAsFactors = FALSE
  )

  # Test process_metadata
  result <- scCulturePredict:::process_metadata(seurat_obj, metadata, barcodes)

  # Check that metadata was added
  expect_s4_class(result, "Seurat")
  expect_true("condition" %in% colnames(result@meta.data))
  expect_true("batch" %in% colnames(result@meta.data))
  expect_equal(result@meta.data$condition, metadata$condition)
  expect_equal(result@meta.data$batch, metadata$batch)
})

test_that("process_metadata handles mismatched barcodes", {
  skip_if_not_installed("Seurat")
  library(Seurat)

  # Create a Seurat object
  counts <- matrix(rpois(100, 5), nrow = 10, ncol = 10)
  colnames(counts) <- paste0("Cell", 1:10)
  rownames(counts) <- paste0("Gene", 1:10)
  seurat_obj <- CreateSeuratObject(counts = counts)

  # Create metadata with some mismatched barcodes
  metadata <- data.frame(
    barcode = paste0("Cell", 1:10),
    value = 1:10,
    stringsAsFactors = FALSE
  )
  rownames(metadata) <- metadata$barcode

  # Create barcodes with some that don't match
  barcodes <- data.frame(
    barcode = c(paste0("Cell", 1:5), paste0("NewCell", 6:10)),
    stringsAsFactors = FALSE
  )

  # Should warn about mismatched barcodes
  expect_warning(
    scCulturePredict:::process_metadata(seurat_obj, metadata, barcodes),
    "barcodes do not match"
  )
})
