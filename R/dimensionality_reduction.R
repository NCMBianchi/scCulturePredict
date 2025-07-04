#' Perform PCA analysis
#'
#' @description
#' Performs Principal Component Analysis (PCA) on a Seurat object and returns
#' the results along with variance explained by each principal component.
#'
#' @param seurat_object A Seurat object containing preprocessed single-cell data.
#' @param n_pcs Integer specifying the number of principal components to compute.
#'   Default is 50.
#' @param features Character vector specifying which features to use for PCA.
#'   If NULL, uses variable features. Default is NULL.
#'
#' @return A list containing:
#' \itemize{
#'   \item seurat_object: The Seurat object with PCA results
#'   \item variance_explained: Data frame with variance explained by each PC
#'   \item elbow_plot: ggplot object showing the elbow plot
#' }
#' @export
#'
#' @examples
#' pca_results <- perform_pca(seurat_object, n_pcs = 50)
perform_pca <- function(seurat_object, n_pcs = 50, features = NULL) {
    # Run PCA
    seurat_object <- Seurat::RunPCA(seurat_object,
        features = features %||% Seurat::VariableFeatures(seurat_object),
        npcs = n_pcs,
        verbose = FALSE
    )

    # Calculate variance explained
    variance_explained <- data.frame(
        PC = seq_len(n_pcs),
        Variance = seurat_object@reductions$pca@stdev^2,
        Cumulative = cumsum(seurat_object@reductions$pca@stdev^2) / sum(seurat_object@reductions$pca@stdev^2)
    )

    # Create elbow plot
    elbow_plot <- ggplot2::ggplot(variance_explained, ggplot2::aes(x = PC, y = Variance)) +
        ggplot2::geom_point() +
        ggplot2::geom_line() +
        ggplot2::theme_minimal() +
        ggplot2::labs(
            title = "PCA Elbow Plot",
            x = "Principal Component",
            y = "Variance Explained"
        )

    return(list(
        seurat_object = seurat_object,
        variance_explained = variance_explained,
        elbow_plot = elbow_plot
    ))
}

#' Run UMAP dimensionality reduction
#'
#' @description
#' Performs Uniform Manifold Approximation and Projection (UMAP) on a Seurat object
#' using PCA results as input. UMAP is a dimensionality reduction technique that
#' preserves both local and global structure of the data.
#'
#' @param seurat_obj Seurat object. The Seurat object to run UMAP on.
#' @param n_neighbors Integer. Number of neighbors for UMAP (default: 30).
#'   This parameter controls how UMAP balances local versus global structure.
#'   Lower values preserve more local structure, while higher values preserve
#'   more global structure.
#' @param min_dist Numeric. Minimum distance for UMAP (default: 0.3).
#'   This parameter controls how tightly UMAP is allowed to pack points together.
#'   Lower values allow points to be closer together, while higher values
#'   force points to be more spread out.
#' @param n_components Integer. Number of UMAP components to compute (default: 2).
#'   This determines the dimensionality of the output embedding.
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A Seurat object with UMAP coordinates added to the reductions slot.
#'   The UMAP coordinates can be accessed using:
#'   \itemize{
#'     \item \code{seurat_obj@reductions$umap} - The UMAP reduction object
#'     \item \code{seurat_obj@reductions$umap@cell.embeddings} - The UMAP coordinates
#'   }
#'
#' @details
#' The function performs the following steps:
#' \enumerate{
#'   \item Checks if PCA has been run, runs it if necessary
#'   \item Performs UMAP using the first 30 principal components
#'   \item Stores the results in the Seurat object
#' }
#'
#' The UMAP algorithm:
#' \itemize{
#'   \item Constructs a high dimensional graph on the data
#'   \item Optimizes a low dimensional layout to preserve the graph structure
#'   \item Can capture both local and global structure in the data
#' }
#'
#' @examples
#' # Example with mock Seurat object
#' library(Seurat)
#' # Create minimal mock data
#' counts <- matrix(rpois(1000, 5), nrow = 100)
#' rownames(counts) <- paste0("Gene", seq_len(100))
#' colnames(counts) <- paste0("Cell", seq_len(10))
#' seurat_obj <- CreateSeuratObject(counts = counts)
#' seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
#' seurat_obj <- FindVariableFeatures(seurat_obj, verbose = FALSE)
#' seurat_obj <- ScaleData(seurat_obj, verbose = FALSE)
#' seurat_obj <- RunPCA(seurat_obj, npcs = 5, verbose = FALSE)
#'
#' # Run UMAP
#' seurat_with_umap <- run_umap(seurat_obj, n_components = 2, n_neighbors = 5)
#' @seealso
#' \code{\link{perform_pca}} for PCA analysis
#' \code{\link{perform_tsne}} for t-SNE analysis
#'
#' @export
run_umap <- function(seurat_obj,
                     n_neighbors = 30,
                     min_dist = 0.3,
                     n_components = 2,
                     verbose = TRUE) {
    # Input validation
    if (!inherits(seurat_obj, "Seurat")) {
        stop("seurat_obj must be a Seurat object")
    }

    if (!is.numeric(n_neighbors) || length(n_neighbors) != 1 || n_neighbors < 0) {
        stop("n_neighbors must be a single positive number")
    }

    if (!is.numeric(min_dist) || length(min_dist) != 1 || min_dist < 0 || min_dist > 1) {
        stop("min_dist must be a single number between 0 and 1")
    }

    if (!is.numeric(n_components) || length(n_components) != 1 || n_components < 1) {
        stop("n_components must be a single positive integer")
    }

    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("verbose must be a single logical value")
    }

    # Check if PCA has been run
    if (!"pca" %in% names(seurat_obj@reductions)) {
        if (verbose) message("PCA not found, running PCA first...")
        seurat_obj <- tryCatch(
            {
                Seurat::RunPCA(seurat_obj,
                    features = Seurat::VariableFeatures(seurat_obj),
                    verbose = verbose
                )
            },
            error = function(e) {
                stop(sprintf("Error running PCA: %s", e$message))
            }
        )
    }

    # Run UMAP
    if (verbose) message("Running UMAP...")
    seurat_obj <- tryCatch(
        {
            Seurat::RunUMAP(seurat_obj,
                dims = 1:30,
                n.neighbors = n_neighbors,
                min.dist = min_dist,
                n.components = n_components,
                verbose = verbose
            )
        },
        error = function(e) {
            stop(sprintf("Error running UMAP: %s", e$message))
        }
    )

    if (verbose) {
        message(sprintf("UMAP complete. Added %d UMAP components", n_components))
    }

    return(seurat_obj)
}

#' Perform t-SNE dimensionality reduction
#'
#' @description
#' Performs t-Distributed Stochastic Neighbor Embedding (t-SNE) on a Seurat object
#' using PCA results as input.
#'
#' @param seurat_object A Seurat object containing PCA results.
#' @param dims Integer vector specifying which PCA dimensions to use.
#'   Default is seq_len(40).
#' @param perplexity Numeric value specifying the perplexity parameter.
#'   Default is 30.
#' @param max_iter Integer specifying the maximum number of iterations.
#'   Default is 1000.
#'
#' @return A Seurat object with t-SNE coordinates added to the metadata.
#' @export
#'
#' @examples
#' seurat_object <- perform_tsne(seurat_object, dims = seq_len(40))
perform_tsne <- function(seurat_object, dims = 1:40, perplexity = 30,
                         max_iter = 1000) {
    # Run t-SNE
    seurat_object <- Seurat::RunTSNE(seurat_object,
        reduction = "pca",
        dims = dims,
        perplexity = perplexity,
        max_iter = max_iter
    )

    # Add t-SNE coordinates to metadata
    tsne_data <- Seurat::FetchData(seurat_object, c("tSNE_1", "tSNE_2"))
    seurat_object@meta.data$tSNE_1 <- tsne_data$tSNE_1
    seurat_object@meta.data$tSNE_2 <- tsne_data$tSNE_2

    return(seurat_object)
}

#' Perform dimensionality reduction pipeline
#'
#' @description
#' Performs a complete dimensionality reduction pipeline including PCA, UMAP,
#' and t-SNE on a Seurat object.
#'
#' @param seurat_object A Seurat object containing preprocessed single-cell data.
#' @param n_pcs Integer specifying the number of principal components to compute.
#'   Default is 50.
#' @param dims Integer vector specifying which PCA dimensions to use for UMAP and t-SNE.
#'   Default is seq_len(40).
#' @param perform_tsne Logical indicating whether to perform t-SNE. Default is TRUE.
#'
#' @return A list containing:
#' \itemize{
#'   \item seurat_object: The Seurat object with all dimensionality reduction results
#'   \item pca_results: Results from PCA analysis
#' }
#' @export
#'
#' @examples
#' results <- perform_dimensionality_reduction(seurat_object)
perform_dimensionality_reduction <- function(seurat_object, n_pcs = 50,
                                             dims = 1:40, perform_tsne = TRUE) {
    # Perform PCA
    pca_results <- perform_pca(seurat_object, n_pcs = n_pcs)
    seurat_object <- pca_results$seurat_object

    # Perform UMAP
    seurat_object <- run_umap(seurat_object)

    # Perform t-SNE if requested
    if (perform_tsne) {
        seurat_object <- perform_tsne(seurat_object, dims = dims)
    }

    return(list(
        seurat_object = seurat_object,
        pca_results = pca_results
    ))
}
