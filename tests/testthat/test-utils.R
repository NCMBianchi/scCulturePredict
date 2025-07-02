# Test utility functions

test_that("check_and_install_packages works correctly", {
  # Test with a package that should be installed
  expect_silent(check_and_install_packages("utils"))
  
  # Test with multiple packages
  expect_silent(check_and_install_packages(c("utils", "stats")))
  
  # Test with invalid repository
  expect_warning(check_and_install_packages("nonexistent_package", 
                                          repos = "http://invalid.repo"))
})

test_that("create_dir_if_not_exists works correctly", {
  # Create a temporary directory for testing
  temp_dir <- tempfile("test_dir")
  on.exit(unlink(temp_dir, recursive = TRUE))
  
  # Test creating a new directory
  expect_silent(create_dir_if_not_exists(temp_dir))
  expect_true(dir.exists(temp_dir))
  
  # Test creating a nested directory
  nested_dir <- file.path(temp_dir, "nested", "dir")
  expect_silent(create_dir_if_not_exists(nested_dir, recursive = TRUE))
  expect_true(dir.exists(nested_dir))
  
  # Test creating an existing directory
  expect_silent(create_dir_if_not_exists(temp_dir))
})

test_that("save_object and load_object work correctly", {
  # Create a temporary file for testing
  temp_file <- tempfile("test_object", fileext = ".rds")
  on.exit(unlink(temp_file))
  
  # Test saving and loading a simple object
  test_obj <- list(a = 1, b = "test")
  expect_silent(save_object(test_obj, temp_file))
  loaded_obj <- load_object(temp_file)
  expect_identical(test_obj, loaded_obj)
  
  # Test saving and loading a complex object
  complex_obj <- list(
    data = matrix(1:9, nrow = 3),
    metadata = data.frame(x = 1:3, y = letters[1:3])
  )
  expect_silent(save_object(complex_obj, temp_file))
  loaded_complex <- load_object(temp_file)
  expect_identical(complex_obj, loaded_complex)
  
  # Test error handling for non-existent file
  expect_error(load_object("nonexistent_file.rds"))
})

test_that("format_number works correctly", {
  # Test basic formatting
  expect_equal(format_number(1.23456, digits = 2), "1.23")
  expect_equal(format_number(1.23456, digits = 4), "1.2346")
  
  # Test vector input
  expect_equal(format_number(c(1.23456, 2.34567), digits = 2), 
               c("1.23", "2.35"))
  
  # Test zero digits
  expect_equal(format_number(1.23456, digits = 0), "1")
})

test_that("calculate_percentage works correctly", {
  # Test basic percentage calculation
  expect_equal(calculate_percentage(50, total = 100), "50.0%")
  expect_equal(calculate_percentage(50, total = 100, digits = 2), "50.00%")
  
  # Test vector input
  expect_equal(calculate_percentage(c(25, 75), total = 100), 
               c("25.0%", "75.0%"))
  
  # Test automatic total calculation
  expect_equal(calculate_percentage(c(25, 75)), 
               c("25.0%", "75.0%"))
})

test_that("is_empty works correctly", {
  # Test various empty objects
  expect_true(is_empty(NULL))
  expect_true(is_empty(c()))
  expect_true(is_empty(data.frame()))
  expect_true(is_empty(NA))
  expect_true(is_empty(c(NA, NA)))
  
  # Test non-empty objects
  expect_false(is_empty(1))
  expect_false(is_empty(c(1, 2)))
  expect_false(is_empty(data.frame(x = 1)))
  expect_false(is_empty(c(NA, 1)))
})

test_that("get_file_extension works correctly", {
  # Test various file extensions
  expect_equal(get_file_extension("test.txt"), "txt")
  expect_equal(get_file_extension("test.csv.gz"), "gz")
  expect_equal(get_file_extension("test"), "")
  expect_equal(get_file_extension("test."), "")
})

test_that("validate_file works correctly", {
  # Create a temporary file for testing
  temp_file <- tempfile("test_file", fileext = ".txt")
  writeLines("test", temp_file)
  on.exit(unlink(temp_file))
  
  # Test valid file
  expect_true(validate_file(temp_file))
  expect_true(validate_file(temp_file, extension = "txt"))
  
  # Test invalid extension
  expect_false(validate_file(temp_file, extension = "csv"))
  
  # Test non-existent file
  expect_false(validate_file("nonexistent_file.txt"))
  
  # Test file without extension
  temp_file_no_ext <- tempfile("test_file")
  writeLines("test", temp_file_no_ext)
  on.exit(unlink(temp_file_no_ext))
  expect_true(validate_file(temp_file_no_ext))
}) 