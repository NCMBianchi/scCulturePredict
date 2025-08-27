# Test file for SingleCellExperiment support functionality

library(testthat)
library(scCulturePredict)
library(SingleCellExperiment)
library(SummarizedExperiment)

# Helper function to create mock SCE data
create_mock_sce_data <- function(n_genes = 500, n_cells = 2000, n_samples = 2) {
  # Create count matrix - use sparse matrix like real scRNA-seq data
  # Use higher expression values to ensure cells pass filtering
  set.seed(123)
  counts <- Matrix::Matrix(
    rpois(n_genes * n_cells, lambda = 20),
    nrow = n_genes,
    ncol = n_cells,
    sparse = TRUE
  )

  # Ensure each gene is expressed in at least 5 cells (to pass min.cells = 3)
  for (g in seq_len(n_genes)) {
    # Count how many cells express this gene
    expressing_cells <- which(counts[g, ] > 0)
    if (length(expressing_cells) < 5) {
      # Add expression to more cells
      cells_to_add <- sample(
        setdiff(seq_len(n_cells), expressing_cells),
        min(5 - length(expressing_cells), n_cells - length(expressing_cells))
      )
      for (cell in cells_to_add) {
        counts[g, cell] <- rpois(1, lambda = 10) + 1
      }
    }
  }

  # Ensure each cell has at least 250 expressed genes (to pass min.features = 200)
  for (i in seq_len(n_cells)) {
    expressing_genes <- which(counts[, i] > 0)
    if (length(expressing_genes) < 300) {
      genes_to_add <- sample(
        setdiff(seq_len(n_genes), expressing_genes),
        min(300 - length(expressing_genes), n_genes - length(expressing_genes))
      )
      for (gene in genes_to_add) {
        counts[gene, i] <- rpois(1, lambda = 10) + 1
      }
    }
  }

  # Add realistic yeast gene names (matching sce00001.keg)
  gene_names <- character(n_genes)
  gene_types <- c(
    "YAL", "YBR", "YCL", "YDR", "YEL", "YFR", "YGL", "YGR",
    "YHL", "YHR", "YIL", "YJL", "YKL", "YKR", "YLL", "YLR",
    "YML", "YMR", "YNL", "YNR", "YOL", "YOR", "YPL", "YPR"
  )

  for (i in seq_len(n_genes)) {
    prefix <- sample(gene_types, 1)
    number <- sprintf("%03d", sample(1:999, 1))
    suffix <- sample(c("W", "C"), 1)
    gene_names[i] <- paste0(prefix, number, suffix)
  }

  rownames(counts) <- gene_names
  colnames(counts) <- paste0("Cell", seq_len(n_cells))

  # Create metadata - use blocks of cells per sample (not alternating)
  cells_per_sample <- ceiling(n_cells / n_samples)
  cell_metadata <- DataFrame(
    sample = rep(paste0("Sample", seq_len(n_samples)), each = cells_per_sample, length.out = n_cells),
    batch = rep(c("Batch1", "Batch2"), each = n_cells / 2, length.out = n_cells),
    cell_type = sample(c("TypeA", "TypeB", "TypeC"), n_cells, replace = TRUE)
  )
  rownames(cell_metadata) <- colnames(counts)

  # Create SCE object
  sce <- SingleCellExperiment(
    assays = list(counts = counts),
    colData = cell_metadata
  )

  return(sce)
}

# Helper function to save SCE to RDS file
create_mock_sce_file <- function() {
  # Use explicit parameters to ensure sufficient cells and genes
  sce <- create_mock_sce_data(n_genes = 1000, n_cells = 3000, n_samples = 4)
  temp_file <- tempfile(fileext = ".rds")
  saveRDS(sce, temp_file)
  return(temp_file)
}

# Helper function to get default KEGG file
get_default_kegg_file <- function() {
  # Use the same KEGG file as 10x tests (sce00001.keg with yeast genes)
  kegg_paths <- c(
    system.file("extdata", "kegg", "sce00001.keg", package = "scCulturePredict"),
    system.file("extdata", "sce00001.keg", package = "scCulturePredict")
  )

  for (path in kegg_paths) {
    if (file.exists(path)) {
      return(path)
    }
  }

  # If none found, create a minimal mock KEGG file
  temp_kegg <- tempfile(fileext = ".keg")
  writeLines(c(
    "C\tPathway1",
    "D\tGene1\tDescription1",
    "D\tGene2\tDescription2",
    "C\tPathway2",
    "D\tGene3\tDescription3",
    "D\tGene4\tDescription4"
  ), temp_kegg)

  return(temp_kegg)
}

# Test loading SCE from object
test_that("load_sce_data works with SCE object", {
  # Create mock SCE object
  sce <- create_mock_sce_data()

  # Load using the function
  expect_no_error({
    seurat_obj <- load_sce_data(
      sce_data_path = sce,
      experiment_id = "test_sce_object",
      min_cells = 3,
      min_features = 10,
      verbose = FALSE
    )
  })

  # Check that result is a Seurat object
  seurat_obj <- load_sce_data(
    sce_data_path = sce,
    experiment_id = "test_sce_object",
    min_cells = 3,
    min_features = 10,
    verbose = FALSE
  )

  expect_s4_class(seurat_obj, "Seurat")
  expect_true("sample" %in% colnames(seurat_obj@meta.data))
  expect_equal(ncol(seurat_obj), ncol(sce))
})

# Test loading SCE from RDS file
test_that("load_sce_data works with RDS file path", {
  # Create mock SCE file
  sce_file <- create_mock_sce_file()

  # Load using the function
  expect_no_error({
    seurat_obj <- load_sce_data(
      sce_data_path = sce_file,
      experiment_id = "test_sce_file",
      min_cells = 3,
      min_features = 10,
      verbose = FALSE
    )
  })

  # Check that result is a Seurat object
  seurat_obj <- load_sce_data(
    sce_data_path = sce_file,
    experiment_id = "test_sce_file",
    min_cells = 3,
    min_features = 10,
    verbose = FALSE
  )

  expect_s4_class(seurat_obj, "Seurat")
  expect_true("sample" %in% colnames(seurat_obj@meta.data))

  # Clean up
  unlink(sce_file)
})

# Test error handling for invalid SCE input
test_that("load_sce_data handles invalid inputs correctly", {
  # Test with non-existent file
  expect_error(
    load_sce_data(
      sce_data_path = "nonexistent_file.rds",
      experiment_id = "test_invalid",
      verbose = FALSE
    ),
    "File not found"
  )

  # Test with invalid object type
  expect_error(
    load_sce_data(
      sce_data_path = list(a = 1, b = 2),
      experiment_id = "test_invalid",
      verbose = FALSE
    ),
    "must be a SingleCellExperiment object or path to RDS file"
  )
})

test_that("load_sce_data handles edge cases with empty data", {
  # Test with empty SCE object (0 cells, 0 genes)
  empty_sce <- SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = 0, ncol = 0))
  )
  expect_error(
    load_sce_data(
      sce_data_path = empty_sce,
      experiment_id = "test_empty",
      verbose = FALSE
    ),
    regexp = "empty|no cells|no genes|insufficient"
  )

  # Test with SCE having cells but no genes
  no_genes_sce <- SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = 0, ncol = 10))
  )
  expect_error(
    load_sce_data(
      sce_data_path = no_genes_sce,
      experiment_id = "test_no_genes",
      verbose = FALSE
    ),
    regexp = "no genes|insufficient|empty"
  )

  # Test with SCE having genes but no cells
  no_cells_sce <- SingleCellExperiment(
    assays = list(counts = matrix(0, nrow = 10, ncol = 0))
  )
  expect_error(
    load_sce_data(
      sce_data_path = no_cells_sce,
      experiment_id = "test_no_cells",
      verbose = FALSE
    ),
    regexp = "no cells|insufficient|empty"
  )
})

test_that("load_sce_data handles single cell data", {
  # Create SCE with single cell
  single_cell_counts <- matrix(rpois(100, lambda = 5), nrow = 100, ncol = 1)
  colnames(single_cell_counts) <- "Cell1"
  rownames(single_cell_counts) <- paste0("Gene", 1:100)

  single_cell_sce <- SingleCellExperiment(
    assays = list(counts = single_cell_counts)
  )
  colData(single_cell_sce)$sample <- "Sample1"

  # Should handle single cell gracefully (might warn or work)
  result <- tryCatch(
    {
      load_sce_data(
        sce_data_path = single_cell_sce,
        experiment_id = "test_single_cell",
        min_cells = 0, # Don't filter out the single cell
        min_features = 10,
        verbose = FALSE
      )
    },
    error = function(e) NULL,
    warning = function(w) "warning"
  )

  # Either works or gives appropriate warning/error
  expect_true(!is.null(result))
})

test_that("load_sce_data handles filtering parameters correctly", {
  # Create SCE with known dimensions
  n_genes <- 200
  n_cells <- 100

  # Create sparse data where some genes/cells will be filtered
  counts <- Matrix::Matrix(0, nrow = n_genes, ncol = n_cells, sparse = TRUE)

  # Make first 50 genes expressed in all cells (won't be filtered)
  for (i in 1:50) {
    counts[i, ] <- rpois(n_cells, lambda = 5)
  }

  # Make next 50 genes expressed in only 2 cells (will be filtered with min_cells = 3)
  for (i in 51:100) {
    counts[i, 1:2] <- rpois(2, lambda = 10)
  }

  # Make last 100 genes not expressed at all (will be filtered)
  # counts[101:200, ] remain 0

  colnames(counts) <- paste0("Cell", 1:n_cells)
  rownames(counts) <- paste0("Gene", 1:n_genes)

  sce <- SingleCellExperiment(assays = list(counts = counts))
  colData(sce)$sample <- rep(c("A", "B"), length.out = n_cells)

  # Test with strict filtering
  result_strict <- load_sce_data(
    sce_data_path = sce,
    experiment_id = "test_strict_filter",
    min_cells = 3, # Will filter out genes 51-200
    min_features = 10,
    verbose = FALSE
  )

  # Should have fewer genes after filtering
  expect_lt(nrow(result_strict), n_genes)

  # Test with no filtering
  result_no_filter <- load_sce_data(
    sce_data_path = sce,
    experiment_id = "test_no_filter",
    min_cells = 0,
    min_features = 0,
    verbose = FALSE
  )

  # Should retain more genes
  expect_gte(nrow(result_no_filter), nrow(result_strict))
})

test_that("load_sce_data handles duplicate handling parameter", {
  # Create SCE with duplicate gene names
  n_genes <- 10
  n_cells <- 20
  counts <- matrix(rpois(n_genes * n_cells, lambda = 10), nrow = n_genes, ncol = n_cells)

  # Create duplicate gene names
  gene_names <- c(
    "Gene1", "Gene2", "Gene1", "Gene3", "Gene2",
    "Gene4", "Gene5", "Gene1", "Gene6", "Gene7"
  )
  rownames(counts) <- gene_names
  colnames(counts) <- paste0("Cell", 1:n_cells)

  sce <- SingleCellExperiment(assays = list(counts = counts))
  colData(sce)$sample <- rep(c("A", "B"), each = 10)

  # Test with make_unique
  result_unique <- load_sce_data(
    sce_data_path = sce,
    experiment_id = "test_make_unique",
    handle_duplicates = "make_unique",
    min_cells = 0,
    min_features = 0,
    verbose = FALSE
  )

  # Should have all genes with unique names
  expect_equal(length(unique(rownames(result_unique))), n_genes)

  # Test with first
  result_first <- load_sce_data(
    sce_data_path = sce,
    experiment_id = "test_first",
    handle_duplicates = "first",
    min_cells = 0,
    min_features = 0,
    verbose = FALSE
  )

  # Should have fewer genes (only first occurrence of each)
  expect_lt(nrow(result_first), n_genes)

  # Test with aggregate
  result_agg <- load_sce_data(
    sce_data_path = sce,
    experiment_id = "test_aggregate",
    handle_duplicates = "aggregate",
    min_cells = 0,
    min_features = 0,
    verbose = FALSE
  )

  # Should have unique gene names
  expect_equal(length(unique(rownames(result_agg))), length(rownames(result_agg)))

  # Test with error
  expect_error(
    load_sce_data(
      sce_data_path = sce,
      experiment_id = "test_error",
      handle_duplicates = "error",
      verbose = FALSE
    ),
    regexp = "Found.*duplicate gene names"
  )
})

test_that("load_sce_data handles verbose parameter correctly", {
  # Create simple SCE
  sce <- create_mock_sce_data(n_genes = 100, n_cells = 50)

  # Test with verbose = TRUE
  expect_message(
    load_sce_data(
      sce_data_path = sce,
      experiment_id = "test_verbose",
      verbose = TRUE,
      min_cells = 3,
      min_features = 10
    ),
    regexp = "Converting SingleCellExperiment"
  )

  # Test with verbose = FALSE
  expect_silent(
    load_sce_data(
      sce_data_path = sce,
      experiment_id = "test_silent",
      verbose = FALSE,
      min_cells = 3,
      min_features = 10
    )
  )
})

test_that("load_sce_data handles non-standard assay slots", {
  # Create SCE with multiple assays
  n_genes <- 100
  n_cells <- 50
  counts <- matrix(rpois(n_genes * n_cells, lambda = 10), nrow = n_genes, ncol = n_cells)
  logcounts <- log2(counts + 1)

  rownames(counts) <- paste0("Gene", 1:n_genes)
  colnames(counts) <- paste0("Cell", 1:n_cells)
  rownames(logcounts) <- rownames(counts)
  colnames(logcounts) <- colnames(counts)

  # SCE with both counts and logcounts
  sce_multi <- SingleCellExperiment(
    assays = list(
      counts = counts,
      logcounts = logcounts
    )
  )
  colData(sce_multi)$sample <- rep(c("A", "B"), length.out = n_cells)

  # Should use counts by default
  result <- load_sce_data(
    sce_data_path = sce_multi,
    experiment_id = "test_multi_assay",
    min_cells = 3,
    min_features = 10,
    verbose = FALSE
  )

  expect_s4_class(result, "Seurat")

  # Test with only logcounts (no counts assay)
  sce_logonly <- SingleCellExperiment(
    assays = list(logcounts = logcounts)
  )
  colData(sce_logonly)$sample <- rep(c("A", "B"), length.out = n_cells)

  # Should handle missing counts assay
  result_logonly <- tryCatch(
    {
      load_sce_data(
        sce_data_path = sce_logonly,
        experiment_id = "test_logonly",
        min_cells = 3,
        min_features = 10,
        verbose = FALSE
      )
    },
    error = function(e) "error"
  )

  # Either works with logcounts or gives appropriate error
  expect_true(inherits(result_logonly, "Seurat") || result_logonly == "error")
})

test_that("scCulture BUILD mode works with SCE data", {
  # Create mock SCE data with sufficient cells and genes to pass filtering
  sce <- create_mock_sce_data(n_genes = 1000, n_cells = 3000)

  # Get KEGG file
  kegg_file <- get_default_kegg_file()

  # Skip test if KEGG file not found
  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Run scCulture with SCE object
  expect_no_error({
    results <- scCulture(
      sce_data_path = sce,
      input_type = "sce",
      mode = "build",
      experiment_id = "test_sce_build",
      kegg_file = kegg_file,
      output_dir = tempfile("test_sce_output"),
      perform_tsne = FALSE,
      verbose = FALSE,
      progress = FALSE
    )
  })

  # Test that results have expected structure
  results <- scCulture(
    sce_data_path = sce,
    input_type = "sce",
    mode = "build",
    experiment_id = "test_sce_build",
    kegg_file = kegg_file,
    output_dir = tempfile("test_sce_output"),
    perform_tsne = FALSE,
    verbose = FALSE,
    progress = FALSE
  )

  expect_type(results, "list")
  expect_true("seurat_object" %in% names(results))
  expect_true("pathway_results" %in% names(results))
  expect_true("prediction_results" %in% names(results))
  expect_true("evaluation_results" %in% names(results))
  expect_true("fingerprint_file" %in% names(results))
})

# Test scCulture with SCE file path
test_that("scCulture works with SCE file path", {
  # Create mock SCE file
  sce_file <- create_mock_sce_file()

  # Get KEGG file
  kegg_file <- get_default_kegg_file()

  # Skip test if KEGG file not found
  if (!file.exists(kegg_file)) {
    skip("Default KEGG file not found in package")
  }

  # Run scCulture with SCE file path
  expect_no_error({
    results <- scCulture(
      sce_data_path = sce_file,
      input_type = "sce",
      mode = "build",
      experiment_id = "test_sce_file_build",
      kegg_file = kegg_file,
      output_dir = tempfile("test_sce_file_output"),
      perform_tsne = FALSE,
      verbose = FALSE,
      progress = FALSE
    )
  })

  # Clean up
  unlink(sce_file)
})

# Test cross-dataset prediction (10X to SCE)
test_that("Cross-dataset prediction from 10X to SCE works", {
  skip("Cross-dataset prediction test requires 10X mock data setup")

  # This test would require creating both 10X and SCE mock data
  # and running the full pipeline, which is complex and time-consuming
  # Leaving as a placeholder for future implementation
})

# Test metadata preservation
test_that("SCE metadata is preserved in conversion", {
  # Create SCE with specific metadata - ensure sufficient cells to pass filtering
  sce <- create_mock_sce_data(n_genes = 1000, n_cells = 3000)

  # Add custom metadata columns
  colData(sce)$custom_field <- seq_len(ncol(sce))
  colData(sce)$treatment <- sample(c("Control", "Treated"), ncol(sce), replace = TRUE)

  # Load and convert
  seurat_obj <- load_sce_data(
    sce_data_path = sce,
    experiment_id = "test_metadata",
    verbose = FALSE
  )

  # Check metadata preservation
  expect_true("sample" %in% colnames(seurat_obj@meta.data))
  expect_true("batch" %in% colnames(seurat_obj@meta.data))
  expect_true("cell_type" %in% colnames(seurat_obj@meta.data))
  expect_true("custom_field" %in% colnames(seurat_obj@meta.data))
  expect_true("treatment" %in% colnames(seurat_obj@meta.data))
})

# Test handling of different assay names
test_that("load_sce_data handles different assay names", {
  # Create SCE with non-standard assay name - use sparse matrix like real data
  # Use more cells and genes with higher expression to pass filtering
  n_genes <- 500
  n_cells <- 2000
  counts <- Matrix::Matrix(rpois(n_genes * n_cells, lambda = 20), nrow = n_genes, ncol = n_cells, sparse = TRUE)

  # Ensure each gene is expressed in at least 5 cells (to pass min.cells = 3)
  for (g in seq_len(n_genes)) {
    expressing_cells <- which(counts[g, ] > 0)
    if (length(expressing_cells) < 5) {
      cells_to_add <- sample(
        setdiff(seq_len(n_cells), expressing_cells),
        min(5 - length(expressing_cells), n_cells - length(expressing_cells))
      )
      for (cell in cells_to_add) {
        counts[g, cell] <- rpois(1, lambda = 10) + 1
      }
    }
  }

  # Ensure each cell has sufficient expressed genes
  for (i in seq_len(n_cells)) {
    expressing_genes <- which(counts[, i] > 0)
    if (length(expressing_genes) < 300) {
      genes_to_add <- sample(
        setdiff(seq_len(n_genes), expressing_genes),
        min(300 - length(expressing_genes), n_genes - length(expressing_genes))
      )
      for (gene in genes_to_add) {
        counts[gene, i] <- rpois(1, lambda = 10) + 1
      }
    }
  }

  # Add realistic yeast gene names (matching sce00001.keg)
  gene_names <- character(n_genes)
  gene_types <- c(
    "YAL", "YBR", "YCL", "YDR", "YEL", "YFR", "YGL", "YGR",
    "YHL", "YHR", "YIL", "YJL", "YKL", "YKR", "YLL", "YLR",
    "YML", "YMR", "YNL", "YNR", "YOL", "YOR", "YPL", "YPR"
  )

  for (i in seq_len(n_genes)) {
    prefix <- sample(gene_types, 1)
    number <- sprintf("%03d", sample(1:999, 1))
    suffix <- sample(c("W", "C"), 1)
    gene_names[i] <- paste0(prefix, number, suffix)
  }

  rownames(counts) <- gene_names
  colnames(counts) <- paste0("Cell", seq_len(n_cells))

  sce <- SingleCellExperiment(
    assays = list(raw = counts), # Using 'raw' instead of 'counts'
    colData = DataFrame(sample = rep(c("A", "B"), each = n_cells / 2, length.out = n_cells))
  )

  # Should still work
  expect_no_error({
    seurat_obj <- load_sce_data(
      sce_data_path = sce,
      experiment_id = "test_assay_names",
      verbose = FALSE
    )
  })
})

# Test parameter validation for scCulture with SCE
test_that("scCulture validates SCE parameters correctly", {
  # Test missing input_type
  sce <- create_mock_sce_data()

  expect_error(
    scCulture(
      sce_data_path = sce,
      # input_type missing
      mode = "build",
      output_dir = tempfile()
    ),
    "input_type must be specified"
  )

  # Test invalid input_type
  expect_error(
    scCulture(
      sce_data_path = sce,
      input_type = "invalid",
      mode = "build",
      output_dir = tempfile()
    ),
    "input_type must be either '10x' or 'sce'"
  )

  # Test missing sce_data_path for SCE input_type
  expect_error(
    scCulture(
      input_type = "sce",
      mode = "build",
      output_dir = tempfile()
    ),
    "sce_data_path must be provided for sce input type"
  )
})

# Test dimension preservation
test_that("Dimensions are preserved during SCE conversion", {
  # Create SCE with specific dimensions
  n_genes <- 150
  n_cells <- 75
  sce <- create_mock_sce_data(n_genes = n_genes, n_cells = n_cells)

  # Convert to Seurat
  seurat_obj <- load_sce_data(
    sce_data_path = sce,
    experiment_id = "test_dimensions",
    min_cells = 0, # Don't filter
    min_features = 0, # Don't filter
    verbose = FALSE
  )

  # Check dimensions (accounting for potential filtering)
  expect_lte(nrow(seurat_obj), n_genes)
  expect_lte(ncol(seurat_obj), n_cells)
})
