# Test visualization functions

# Create a mock Seurat object for testing
create_mock_seurat <- function() {
  # Create a small Seurat object with UMAP coordinates
  counts <- matrix(rpois(5000, 10), nrow = 100, ncol = 50)
  rownames(counts) <- paste0("gene", 1:100)
  colnames(counts) <- paste0("cell", 1:50)

  seurat <- Seurat::CreateSeuratObject(counts = counts)

  # Add UMAP coordinates
  seurat@meta.data$UMAP_1 <- rnorm(50)
  seurat@meta.data$UMAP_2 <- rnorm(50)
  seurat@meta.data$sample <- rep(c("A", "B"), each = 25)

  return(seurat)
}

# Create mock evaluation results
create_mock_evaluation <- function() {
  list(
    direct_accuracy = data.frame(
      sample = c("A", "B"),
      correct = c(4, 3),
      percent = c("80%", "60%")
    ),
    direct_table = table(
      Actual = c("A", "A", "B", "B"),
      Predicted = c("A", "B", "A", "B")
    ),
    svm_accuracy = data.frame(
      sample = c("A", "B"),
      correct = c(5, 4),
      percent = c("100%", "80%")
    ),
    svm_table = table(
      Actual = c("A", "A", "B", "B"),
      Predicted = c("A", "B", "A", "B")
    )
  )
}

test_that("create_umap_plot works correctly", {
  seurat <- create_mock_seurat()

  # Test basic plot creation
  p <- create_umap_plot(seurat)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "UMAP Visualization")

  # Test custom parameters
  p <- create_umap_plot(seurat,
    color_by = "sample",
    point_size = 1.5,
    title = "Custom Title",
    legend_title = "Custom Legend"
  )
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "Custom Title")
  expect_equal(p$labels$colour, "Custom Legend")
})

test_that("create_accuracy_plot works correctly", {
  eval_results <- create_mock_evaluation()

  # Test direct method plot
  p <- create_accuracy_plot(eval_results, method = "direct")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "DIRECT Prediction Accuracy")

  # Test SVM method plot
  p <- create_accuracy_plot(eval_results, method = "svm")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "SVM Prediction Accuracy")

  # Test custom title
  p <- create_accuracy_plot(eval_results, method = "direct", title = "Custom Title")
  expect_equal(p$labels$title, "Custom Title")
})

test_that("create_confusion_heatmap works correctly", {
  eval_results <- create_mock_evaluation()

  # Test direct method heatmap
  p <- create_confusion_heatmap(eval_results, method = "direct")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "DIRECT Confusion Matrix")

  # Test SVM method heatmap
  p <- create_confusion_heatmap(eval_results, method = "svm")
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$title, "SVM Confusion Matrix")

  # Test custom title
  p <- create_confusion_heatmap(eval_results, method = "direct", title = "Custom Title")
  expect_equal(p$labels$title, "Custom Title")
})

test_that("save_visualization_plots works correctly", {
  seurat <- create_mock_seurat()
  eval_results <- create_mock_evaluation()

  # Create temporary directory for testing
  temp_dir <- tempfile("test_plots")
  dir.create(temp_dir)
  on.exit(unlink(temp_dir, recursive = TRUE))

  # Test saving plots
  expect_silent(save_visualization_plots(seurat, eval_results, temp_dir, verbose = FALSE))

  # Check if files were created
  expected_files <- c(
    "plot_umap.png",
    "plot_direct_accuracy.png",
    "plot_svm_accuracy.png",
    "plot_direct_confusion.png",
    "plot_svm_confusion.png"
  )

  for (file in expected_files) {
    expect_true(file.exists(file.path(temp_dir, file)))
  }

  # Test with custom prefix
  expect_silent(save_visualization_plots(seurat, eval_results, temp_dir, prefix = "custom", verbose = FALSE))
  expect_true(file.exists(file.path(temp_dir, "custom_umap.png")))
})
