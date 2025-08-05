# Test dimensionality reduction functions

# Create a mock Seurat object for testing
create_mock_seurat <- function() {
  # Create a small Seurat object with variable features
  counts <- matrix(rpois(5000, 10), nrow = 100, ncol = 50)
  rownames(counts) <- paste0("gene", 1:100)
  colnames(counts) <- paste0("cell", 1:50)

  seurat <- Seurat::CreateSeuratObject(counts = counts)
  seurat <- Seurat::NormalizeData(seurat)
  seurat <- Seurat::FindVariableFeatures(seurat)
  seurat <- Seurat::ScaleData(seurat)

  return(seurat)
}

test_that("perform_pca works correctly", {
  seurat <- create_mock_seurat()

  # Test basic PCA
  results <- perform_pca(seurat, n_pcs = 10)
  expect_s4_class(results$seurat_object, "Seurat")
  expect_s3_class(results$variance_explained, "data.frame")
  expect_s3_class(results$elbow_plot, "ggplot")

  # Check variance explained data frame
  expect_equal(nrow(results$variance_explained), 10)
  expect_true(all(c("PC", "Variance", "Cumulative") %in%
    names(results$variance_explained)))

  # Test with custom features
  custom_features <- paste0("gene", 1:20)
  results <- perform_pca(seurat, n_pcs = 5, features = custom_features)
  expect_equal(nrow(results$variance_explained), 5)
})

test_that("perform_umap works correctly", {
  seurat <- create_mock_seurat()

  # Run PCA first
  seurat <- Seurat::RunPCA(seurat, features = Seurat::VariableFeatures(seurat))

  # Test basic UMAP
  result <- perform_umap(seurat, dims = 1:5)
  expect_s4_class(result, "Seurat")
  expect_true(all(c("UMAP_1", "UMAP_2") %in% names(result@meta.data)))

  # Test with custom parameters
  result <- perform_umap(seurat,
    dims = 1:3,
    n_neighbors = 5,
    min_dist = 0.1,
    metric = "euclidean"
  )
  expect_s4_class(result, "Seurat")
})

test_that("perform_tsne works correctly", {
  seurat <- create_mock_seurat()

  # Run PCA first
  seurat <- Seurat::RunPCA(seurat, features = Seurat::VariableFeatures(seurat))

  # Test basic t-SNE
  result <- perform_tsne(seurat, dims = 1:5)
  expect_s4_class(result, "Seurat")
  expect_true(all(c("tSNE_1", "tSNE_2") %in% names(result@meta.data)))

  # Test with custom parameters
  result <- perform_tsne(seurat,
    dims = 1:3,
    perplexity = 5,
    max_iter = 500
  )
  expect_s4_class(result, "Seurat")
})

test_that("perform_dimensionality_reduction works correctly", {
  seurat <- create_mock_seurat()

  # Test complete pipeline without t-SNE
  results <- perform_dimensionality_reduction(seurat,
    n_pcs = 10,
    dims = 1:5,
    perform_tsne = FALSE
  )
  expect_s4_class(results$seurat_object, "Seurat")
  expect_s3_class(results$pca_results, "list")
  expect_true(all(c("UMAP_1", "UMAP_2") %in%
    names(results$seurat_object@meta.data)))
  expect_false(all(c("tSNE_1", "tSNE_2") %in%
    names(results$seurat_object@meta.data)))

  # Test complete pipeline with t-SNE
  results <- perform_dimensionality_reduction(seurat,
    n_pcs = 10,
    dims = 1:5,
    perform_tsne = TRUE
  )
  expect_s4_class(results$seurat_object, "Seurat")
  expect_true(all(c("UMAP_1", "UMAP_2", "tSNE_1", "tSNE_2") %in%
    names(results$seurat_object@meta.data)))
})
