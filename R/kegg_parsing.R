#' Get best available data layer from Seurat object
#'
#' @description
#' Helper function to determine the best available data layer for FetchData.
#' Checks for normalized data layers and falls back to counts if needed.
#' Handles both Seurat v4 and v5 layer structures.
#'
#' @param seurat_object A Seurat object
#' @return Character string specifying the best layer to use
#' @keywords internal
#' @examples
#' # Example with mock KEG file content
#' keg_content <- c(
#'     "ENTRY       hsa00010                    Pathway",
#'     "NAME        Glycolysis / Gluconeogenesis",
#'     "GENE        1234  HK1; hexokinase 1",
#'     "GENE        5678  GPI; glucose-6-phosphate isomerase",
#'     "///"
#' )
#' temp_file <- tempfile(fileext = ".keg")
#' writeLines(keg_content, temp_file)
#' kegg_data <- parse_kegg_keg(temp_file)
#' print(kegg_data)
#' unlink(temp_file)
get_best_data_layer <- function(seurat_object) {
    # Get the default assay (usually "RNA")
    default_assay <- Seurat::DefaultAssay(seurat_object)

    # Check if the assay exists
    if (!default_assay %in% names(seurat_object@assays)) {
        stop(sprintf("Default assay '%s' not found in Seurat object", default_assay))
    }

    # Try to access normalized data layer
    tryCatch(
        {
            data_layer <- seurat_object[[default_assay]]$data
            if (!is.null(data_layer) && nrow(data_layer) > 0) {
                return("data")
            }
        },
        warning = function(w) {
            # Suppress layer warnings and continue to next option
        },
        error = function(e) {
            # Continue to next option if data layer access fails
        }
    )

    # Try to access counts layer as fallback
    tryCatch(
        {
            counts_layer <- seurat_object[[default_assay]]$counts
            if (!is.null(counts_layer) && nrow(counts_layer) > 0) {
                return("counts")
            }
        },
        warning = function(w) {
            # Suppress layer warnings and continue to next option
        },
        error = function(e) {
            # Continue to next option if counts layer access fails
        }
    )

    # Final fallback - let Seurat decide (will likely use default layer)
    return(NULL)
}

#' Parse KEGG pathway file
#'
#' @description
#' Parses a KEGG pathway file (.keg) and extracts gene names for each pathway.
#'
#' @param file_path Character string specifying the path to the KEGG pathway file.
#' Default is "sce00001.keg".
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A list where each element is a character vector of gene names for a pathway.
#' @export
#'
#' @examples
#' # Example with mock KEG file
#' # Create a temporary KEG file
#' keg_content <- c(
#'     "ENTRY       hsa00010                    Pathway",
#'     "NAME        Glycolysis / Gluconeogenesis - Homo sapiens (human)",
#'     "DESCRIPTION Glycolysis is the process of converting glucose into pyruvate",
#'     "CLASS       Metabolism; Carbohydrate metabolism",
#'     "PATHWAY_MAP hsa00010  Glycolysis / Gluconeogenesis",
#'     "GENE        5213  PFKM; phosphofructokinase, muscle [KO:K00850] [EC:2.7.1.11]",
#'     "GENE        5214  PFKL; phosphofructokinase, liver [KO:K00850] [EC:2.7.1.11]",
#'     "GENE        5211  PFKP; phosphofructokinase, platelet [KO:K00850] [EC:2.7.1.11]",
#'     "///"
#' )
#'
#' temp_file <- tempfile(fileext = ".keg")
#' writeLines(keg_content, temp_file)
#'
#' # Parse the KEG file
#' kegg_data <- parse_kegg_keg(temp_file)
#' print(names(kegg_data))
#' print(kegg_data$genes[1:3])
#'
#' # Clean up
#' unlink(temp_file)
parse_kegg_keg <- function(file_path = "sce00001.keg", verbose = TRUE) {
    lines <- readLines(file_path)
    pathways <- list()
    current_pathway <- NULL

    for (line in lines) {
        if (startsWith(line, "C")) {
            current_pathway <- sub("^C\\s+\\d+\\s+", "", line)
            pathways[[current_pathway]] <- c()
        } else if (startsWith(line, "D")) {
            matches <- regmatches(line, gregexpr("Y[A-Z][LR][0-9]+[WC]", line))
            if (length(matches[[1]]) > 0) {
                gene_names <- matches[[1]]
                pathways[[current_pathway]] <- c(pathways[[current_pathway]], gene_names)
            }
        }
    }

    return(pathways)
}

#' Build transcriptomic fingerprints
#'
#' @description
#' Builds transcriptomic fingerprints by aggregating gene expression by pathway
#' and creating signature profiles for each condition.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param kegg_pathways A list of KEGG pathways and their associated genes.
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#' \itemize{
#'   \item pathway_matrix: Matrix of pathway expression values
#'   \item signature_matrix: Matrix of condition-specific signature profiles
#' }
#' @export
#'
#' @examples
#' # Example with minimal mock data
#' library(Seurat)
#' # Create small mock dataset
#' counts <- matrix(rpois(2000, 5), nrow = 200)
#' rownames(counts) <- paste0("Gene", seq_len(200))
#' colnames(counts) <- paste0("Cell", seq_len(10))
#' metadata <- data.frame(
#'     row.names = colnames(counts),
#'     culture = rep(c("TypeA", "TypeB"), each = 5)
#' )
#' seurat_obj <- CreateSeuratObject(counts = counts, meta.data = metadata)
#'
#' # Create mock pathway database
#' pathways <- list(
#'     Pathway1 = sample(rownames(counts), 10),
#'     Pathway2 = sample(rownames(counts), 15),
#'     Pathway3 = sample(rownames(counts), 12)
#' )
#'
#' # Build fingerprints
#' fingerprints <- build_fingerprints(
#'     seurat_obj,
#'     group_by = "culture",
#'     pathways = pathways
#' )
build_fingerprints <- function(seurat_object, kegg_pathways, verbose = TRUE) {
    # Aggregate gene expression by pathway
    pathway_names <- names(kegg_pathways)
    if (verbose) {
        message(sprintf("Processing %d pathways for fingerprint building...", length(pathway_names)))
    }

    # Determine best available data layer
    best_layer <- get_best_data_layer(seurat_object)

    pathway_expression <- lapply(pathway_names, function(pathway) {
        pathway_genes <- kegg_pathways[[pathway]]
        valid_genes <- pathway_genes[!is.na(pathway_genes) & pathway_genes %in% rownames(seurat_object)]
        if (length(valid_genes) == 0) {
            return(rep(NA, ncol(seurat_object)))
        }

        # Use layer parameter only if best_layer is not NULL
        if (is.null(best_layer)) {
            pathway_data <- Seurat::FetchData(seurat_object, vars = valid_genes)
        } else {
            pathway_data <- Seurat::FetchData(seurat_object, vars = valid_genes, layer = best_layer)
        }
        return(rowMeans(pathway_data, na.rm = TRUE))
    })

    # Create matrix safely with proper column names
    pathway_expression_matrix <- do.call(cbind, pathway_expression)
    colnames(pathway_expression_matrix) <- pathway_names

    # Create signature profiles for each condition
    cell_conditions <- seurat_object@meta.data$sample
    unique_conditions <- unique(cell_conditions)
    signature_profiles <- lapply(unique_conditions, function(condition) {
        condition_rows <- which(cell_conditions == condition)
        colMeans(pathway_expression_matrix[condition_rows, ], na.rm = TRUE)
    })
    signature_matrix <- do.call(cbind, signature_profiles)
    colnames(signature_matrix) <- unique_conditions

    return(list(pathway_matrix = pathway_expression_matrix, signature_matrix = signature_matrix))
}

#' Calculate pathway activities only (for PREDICT mode)
#'
#' @description
#' Calculates pathway activity scores for new data without building signature profiles.
#' This function is used in PREDICT mode where signature matrices are already available
#' from pre-built fingerprints.
#'
#' @details
#' This function performs only the pathway activity calculation part of fingerprint
#' building. It aggregates gene expression by pathway but does not create signature
#' profiles for different conditions. This is the appropriate function to use when
#' applying pre-built fingerprints to new unlabeled data.
#'
#' @param seurat_object A Seurat object containing single-cell data.
#' @param kegg_pathways A list of KEGG pathways with gene names.
#' @param verbose Logical. Whether to print progress messages (default: TRUE).
#'
#' @return A matrix of pathway activity scores with cells as rows and pathways as columns.
#' @export
#'
#' @examples
#' # Example with mock data
#' library(Seurat)
#' # Create minimal mock Seurat object
#' counts <- matrix(rpois(1000, 5), nrow = 100)
#' rownames(counts) <- paste0("Gene", seq_len(100))
#' colnames(counts) <- paste0("Cell", seq_len(10))
#' seurat_obj <- CreateSeuratObject(counts = counts)
#' seurat_obj <- NormalizeData(seurat_obj, verbose = FALSE)
#'
#' # Create mock pathways
#' mock_pathways <- list(
#'     Pathway1 = paste0("Gene", seq_len(20)),
#'     Pathway2 = paste0("Gene", 21:40),
#'     Pathway3 = paste0("Gene", 41:60)
#' )
#'
#' # Calculate pathway activities
#' pathway_activities <- calculate_pathway_activities(
#'     seurat_obj, mock_pathways,
#'     verbose = FALSE
#' )
calculate_pathway_activities <- function(seurat_object, kegg_pathways, verbose = TRUE) {
    # Aggregate gene expression by pathway
    pathway_names <- names(kegg_pathways)
    if (verbose) {
        message(sprintf("Calculating pathway activities for %d pathways...", length(pathway_names)))
    }

    # Determine best available data layer
    best_layer <- get_best_data_layer(seurat_object)

    pathway_expression <- lapply(pathway_names, function(pathway) {
        pathway_genes <- kegg_pathways[[pathway]]
        valid_genes <- pathway_genes[!is.na(pathway_genes) & pathway_genes %in% rownames(seurat_object)]
        if (length(valid_genes) == 0) {
            return(rep(NA, ncol(seurat_object)))
        }

        # Use layer parameter only if best_layer is not NULL
        if (is.null(best_layer)) {
            pathway_data <- Seurat::FetchData(seurat_object, vars = valid_genes)
        } else {
            pathway_data <- Seurat::FetchData(seurat_object, vars = valid_genes, layer = best_layer)
        }
        return(rowMeans(pathway_data, na.rm = TRUE))
    })

    # Create matrix safely with proper column names
    pathway_expression_matrix <- do.call(cbind, pathway_expression)
    colnames(pathway_expression_matrix) <- pathway_names

    return(pathway_expression_matrix)
}
