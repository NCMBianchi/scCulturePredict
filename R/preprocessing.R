#' Preprocess single-cell data
#'
#' @description
#' Performs standard preprocessing steps on a Seurat object including normalization,
#' variable feature selection, and scaling. This function implements the standard
#' Seurat workflow for single-cell RNA-seq data preprocessing.
#'
#' @param seurat_obj Seurat object. The Seurat object to preprocess.
#' @param n_features Integer. Number of variable features to select (default: 2000).
#'   These are the most informative genes that will be used for downstream analysis.
#' @param scale_factor Numeric. Scale factor for data scaling (default: 10000).
#'   This is used in the normalization step to scale the data to a common total.
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A processed Seurat object containing:
#' \itemize{
#'   \item Normalized data in the \code{data} slot
#'   \item Scaled data in the \code{scale.data} slot
#'   \item Variable features identified and stored in the object
#' }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Normalizes the data using LogNormalize method
#'   \item Identifies variable features using the vst method
#'   \item Scales the data for all genes
#' }
#'
#' The normalization step:
#' \itemize{
#'   \item Normalizes the gene expression values for each cell
#'   \item Uses a scale factor of 10,000 by default
#'   \item Applies a log transformation
#' }
#'
#' The variable feature selection:
#' \itemize{
#'   \item Uses the vst method to identify highly variable genes
#'   \item Selects the top 2000 genes by default
#'   \item These genes are used for downstream analysis
#' }
#'
#' The scaling step:
#' \itemize{
#'   \item Centers and scales the data
#'   \item Makes genes comparable across cells
#'   \item Prepares data for dimensionality reduction
#' }
#'
#' @examples
#' # Example with mock Seurat object
#' library(Seurat)
#' # Create minimal mock data
#' counts <- matrix(rpois(1000, 5), nrow = 100)
#' rownames(counts) <- paste0("Gene", seq_len(100))
#' colnames(counts) <- paste0("Cell", seq_len(10))
#' metadata <- data.frame(
#'   row.names = colnames(counts),
#'   condition = rep(c("A", "B"), each = 5)
#' )
#' seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
#'
#' # Preprocess the data
#' processed_obj <- preprocess_data(seurat_obj)
#' @seealso
#' \code{\link{load_10x_data}} for loading the data
#' \code{\link{reduce_dimensions}} for dimensionality reduction
#'
#' @export
preprocess_data <- function(seurat_obj,
                            n_features = 2000,
                            scale_factor = 10000,
                            verbose = TRUE) {
  # Input validation
  if (!inherits(seurat_obj, "Seurat")) {
    stop("seurat_obj must be a Seurat object")
  }

  if (!is.numeric(n_features) || length(n_features) != 1 || n_features < 0) {
    stop("n_features must be a single positive number")
  }

  if (!is.numeric(scale_factor) || length(scale_factor) != 1 || scale_factor <= 0) {
    stop("scale_factor must be a single positive number")
  }

  if (!is.logical(verbose) || length(verbose) != 1) {
    stop("verbose must be a single logical value")
  }

  # Check if data is already normalized (avoiding layer warning)
  data_normalized <- tryCatch(
    {
      data_layer <- seurat_obj[["RNA"]]$data
      !is.null(data_layer) && nrow(data_layer) > 0
    },
    warning = function(w) {
      # Suppress layer warnings and assume normalization is needed
      FALSE
    },
    error = function(e) {
      # If error accessing layer, assume normalization is needed
      FALSE
    }
  )

  if (data_normalized) {
    if (verbose) message("Data is already normalized, skipping normalization step")
  } else {
    if (verbose) message("Normalizing data...")
    seurat_obj <- tryCatch(
      {
        Seurat::NormalizeData(seurat_obj,
          normalization.method = "LogNormalize",
          scale.factor = scale_factor
        )
      },
      error = function(e) {
        stop(sprintf("Error normalizing data: %s", e$message))
      }
    )
  }

  # Find variable features
  if (verbose) message("Finding variable features...")
  seurat_obj <- tryCatch(
    {
      Seurat::FindVariableFeatures(seurat_obj,
        selection.method = "vst",
        nfeatures = n_features
      )
    },
    error = function(e) {
      stop(sprintf("Error finding variable features: %s", e$message))
    }
  )

  # Scale data (following original script approach - only variable features)
  if (verbose) message("Scaling data...")
  seurat_obj <- tryCatch(
    {
      Seurat::ScaleData(seurat_obj,
        features = Seurat::VariableFeatures(object = seurat_obj)
      )
    },
    error = function(e) {
      stop(sprintf("Error scaling data: %s", e$message))
    }
  )

  if (verbose) {
    message(sprintf(
      "Preprocessing complete. Selected %d variable features",
      length(Seurat::VariableFeatures(seurat_obj))
    ))
  }

  return(seurat_obj)
}

#' Reduce dimensions of single-cell data
#'
#' @description
#' Performs dimensionality reduction using PCA and UMAP on a Seurat object.
#' Optionally performs t-SNE analysis.
#'
#' @param seurat_object A Seurat object containing preprocessed single-cell data.
#' @param perform_tsne Logical indicating whether to perform t-SNE. Default is TRUE.
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A Seurat object with PCA and UMAP coordinates added to metadata.
#' If perform_tsne is TRUE, t-SNE coordinates are also added.
#' @export
#'
#' @examples
#' # Create mock Seurat object with enough cells for PCA
#' library(Seurat)
#' counts <- matrix(rpois(100000, 5), nrow = 200, ncol = 500)
#' rownames(counts) <- paste0("Gene", seq_len(200))
#' colnames(counts) <- paste0("Cell", seq_len(500))
#' seurat_obj <- CreateSeuratObject(counts = counts)
#'
#' # Preprocess the data first (required for PCA)
#' seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
#' seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
#' seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
#'
#' # Run dimensionality reduction
#' seurat_obj_reduced <- reduce_dimensions(seurat_obj, perform_tsne = FALSE, verbose = FALSE)
#'
#' # Check that UMAP coordinates were added
#' head(seurat_obj_reduced@meta.data[, c("UMAP_1", "UMAP_2")])
reduce_dimensions <- function(seurat_object, perform_tsne = TRUE, verbose = TRUE) {
  if (verbose) message("Running PCA...")
  # Run PCA
  seurat_object <- Seurat::RunPCA(seurat_object, features = Seurat::VariableFeatures(object = seurat_object), verbose = FALSE)

  if (verbose) message("Running UMAP...")
  # Run UMAP with specific parameters from the original study
  seurat_object <- Seurat::RunUMAP(seurat_object, reduction = "pca", dims = 1:40, n.neighbors = 4, min.dist = .05)

  # Add UMAP coordinates to metadata (Seurat v5 compatible)
  umap_embeddings <- seurat_object[["umap"]]@cell.embeddings
  seurat_object@meta.data$UMAP_1 <- umap_embeddings[, 1]
  seurat_object@meta.data$UMAP_2 <- umap_embeddings[, 2]

  # Run t-SNE if requested
  if (perform_tsne) {
    if (verbose) message("Running t-SNE...")
    seurat_object <- Seurat::RunTSNE(seurat_object, reduction = "pca", dims = 1:40)

    # Add t-SNE coordinates to metadata (Seurat v5 compatible)
    tsne_embeddings <- seurat_object[["tsne"]]@cell.embeddings
    seurat_object@meta.data$tSNE_1 <- tsne_embeddings[, 1]
    seurat_object@meta.data$tSNE_2 <- tsne_embeddings[, 2]
  }

  if (verbose) message("Dimensionality reduction complete")

  return(seurat_object)
}
