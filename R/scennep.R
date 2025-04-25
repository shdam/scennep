#' Single-Cell Nearest-Neighbor Pseudobulking of scRNA-seq expression data
#'
#' This function takes a Seurat or SCE object and computes pseudobulk
#' expression profiles for each cell based on its neighbors. The pseudobulk profiles are
#' created by aggregating the expression values of a cell and its top nearest neighbors, specified
#' by the `nn_count`. If a cell has fewer neighbors than `nn_count`, all available neighbors
#' are used.
#'
#' @param obj A Seurat or SCE object or matrix containing single-cell expression data.
#' @param nn_count The number of nearest neighbors to consider for pseudobulking. If a
#'   cell has fewer than this number of neighbors, all are used. Default is 100.
#' @param nn_top The number of top nearest neighbors to use for aggregation.
#' @param nn_cutoff The strength cutoff for shared neighbors.
#' @param FUN Aggregation function. Defaults to rowSums.
#' @param npcs Number of PCs to use for the neighborhood graph.
#' @param markers Subset to markers of interest.
#' @param as If you provide a matrix, decide if you want to use Seurat or Bioconductor methods. Defaults to Seurat.
#' @param assay The name of the assay in the SingleCellExperiment object.
#' @param mc.cores The number of cores to use for parallel processing. By default, it uses
#'   one less than the total number of cores available to prevent the system from locking up.
#' @param pb Logical. If TRUE, a progress bar is printed.
#' @param return_S4 Logical. If TRUE, the Seurat/SCE object is returned with a 'scennep' assay.
#'
#' @return The input Seurat/SCE object with an additional assay named "scennep", which contains
#'   the pseudobulked expression data.
#'   
#' @importFrom parallel mcapply
#' @importFrom Matrix rowSums rowMeans
#'
#' @examples
#' # Assuming `pbmc` is a Seurat object with SCTransform normalization and an SNN graph computed
#' pbmc <- scennep(pbmc, nn_count = 100, mc.cores = 8, return_S4 = TRUE)
#'
#' @export
scennep <- function(
        obj,
        nn_count = 50,
        nn_top = nn_count,
        nn_cutoff = 1/5,
        pc_explained = .90,
        npcs = NULL,
        FUN = NULL,
        markers = NULL,
        as = c("seurat", "bioc"),
        flavor = c("lognormal", "SCT", "CLR", "none"),
        return_S4 = ifelse(inherits(obj, "matrix"), FALSE, TRUE),
        normalize_data = TRUE,
        assay = c("counts", "exprs", "logcounts"),
        pb = FALSE,
        mc.cores = 1) {


    # Stop if
    stopifnot(
        "Valid object types are 'Seurat', 'matrix', 'SingleCellExperiment'"
        = class(obj)[1] %in%
        c("Seurat", "matrix", "dgCMatrix", "SingleCellExperiment"))

    # Match args
    as <- match.arg(as)
    flavor <- match.arg(flavor)
    assay <- match.arg(assay)
  
    if (inherits(obj, "Seurat")) {
        as <- "seurat"
    } else if (inherits(obj, "SummarizedExperiment")) {
        as <- "bioc"
    }
  
    if (as == "seurat") {
        check_package("Seurat")
    } else {
        invisible(sapply(
            c("SingleCellExperiment", "scater", "scuttle", "scran", "igraph"), 
            check_package, repo = "bioc"))
    }

    # Define APPLY
    APPLY <- set_apply(mc.cores, pb)

    # Create object type
    if (class(obj)[1] %in% c("matrix", "dgCMatrix")) {
        if (as == "seurat") {
            obj <- Seurat::CreateSeuratObject(
                counts = Seurat::CreateAssayObject(counts = obj, key = "RNA_"), 
                assay = "RNA")
        } else {
            obj <- SingleCellExperiment::SingleCellExperiment(obj)
            SummarizedExperiment::assayNames(obj) <- assay
        }
    }

    # Normalize data
    if (as == "bioc") {
        if (assay == "counts" & normalize_data) {
            obj <- scuttle::logNormCounts(obj)
            assay <- "logcounts"
        }
        # Build PCA
        obj <- scater::runPCA(obj, scale = TRUE, exprs_values = assay)
        cumulative_variance <- SingleCellExperiment::reducedDim(obj, "PCA") |> 
            attr("percentVar") |> 
            cumsum() / 100
    } else if (as == "seurat") {
        if ("scennep" %in% Seurat::Assays(obj)){
            warning("Assay 'scennep' already exists! It will be overwritten.\n")
            Seurat::DefaultAssay(obj) <- Seurat::Assays(obj)[Seurat::Assays(obj) != "scennep"][1]
            obj[["scennep"]] <- NULL
        }
        # Normalize flavor
        if (flavor == "SCT") {
            if (!"SCT" %in% Seurat::Assays(obj)) {
                obj <- Seurat::SCTransform(obj)
            } else {
                Seurat::DefaultAssay(obj) <- "SCT"
            }
        } else if (flavor == "lognormal") {
            obj <- Seurat::NormalizeData(obj)
        } else if (flavor == "CLR") {
            obj <- Seurat::NormalizeData(obj, normalization.method = "CLR", margin = 2)
        } else if (flavor == "none") {
            obj <- SeuratObject::SetAssayData(obj, layer = "data", new.data = SeuratObject::GetAssayData(obj, layer = "counts"))
        }
        # Build PCA
        if (!"scale.data" %in% SeuratObject::Layers(obj)) obj <- Seurat::ScaleData(obj)
        obj <- Seurat::FindVariableFeatures(obj)
        if (!"pca" %in% Seurat::Reductions(obj)) obj <- Seurat::RunPCA(obj)
        variance <- obj[['pca']]@stdev**2
        cumulative_variance <- cumsum(variance) / sum(variance)
    }
    if (is(npcs, "NULL")) {
        npcs <- which(cumulative_variance >= pc_explained)[1]
        npcs <- ifelse(is.na(npcs), length(cumulative_variance), npcs)
    }

    # Build SNN graph
    snn_graph <- build_snn(obj, nn_count, npcs)
    
    # Extract expression data
    expr <- extract_expression(obj, assay, normalize_data)

    if (is.null(markers)) markers <- rownames(expr)
    expr <- expr[markers, ]
    cols <- colnames(obj)


    # More memory efficient to rm S4 object
    if (!return_S4) {
        rm(obj)
        gc()
    }
    num_cells <- ncol(snn_graph)
    message("Pseudo-bulking each cell with its ", nn_top, " nearest neighbors")
    pseudobulked_expr <- APPLY(seq_len(num_cells), function(cell_id) {
        neighbors <- which(snn_graph[, cell_id] > nn_cutoff)
        nn_top <- min(nn_top, length(neighbors))
        # if (length(neighbors) > nn_top) {
            # Order neighbors by SNN strength
        neighbor_weights <- snn_graph[neighbors, cell_id]
        top_neighbors <- order(neighbor_weights, decreasing = TRUE)[seq_len(nn_top)]
        neighbors <- neighbors[top_neighbors]
        neighbor_weights <- neighbor_weights[top_neighbors]
        # }
        if (length(neighbors) == 1) {
            return(expr[markers, neighbors])
            # neighbors <- c(neighbors, neighbors)
            # neighbor_weights <- c(neighbor_weights, neighbor_weights)
        } 
        if (is.null(FUN)) {
            pseudobulked <- Matrix::rowSums(expr[markers, neighbors] * neighbor_weights) / sum(neighbor_weights)
        } else {
            pseudobulked <- FUN(expr[markers, neighbors])
        }
        return(pseudobulked)
    })

    rm(expr)
    pseudobulked_expr <- do.call(cbind, pseudobulked_expr)
    colnames(pseudobulked_expr) <- cols

    if (!return_S4) return(pseudobulked_expr)

    message("Adding pseudo-bulked expression data to assay 'scennep'")
    if (as == "seurat") {
        if (normalize_data) {
            pseudobulked_expr <- Seurat::CreateAssayObject(data = pseudobulked_expr, key = "scennep_")
        } else {
            pseudobulked_expr <- Seurat::CreateAssayObject(counts = pseudobulked_expr, key = "scennep_")
        }
        
        obj[["scennep"]] <- pseudobulked_expr
        Seurat::VariableFeatures(obj, assay = "scennep") <- Seurat::VariableFeatures(obj, assay = Seurat::DefaultAssay(obj))
        Seurat::DefaultAssay(obj) <- "scennep"
    } else if (as == "bioc"){
        # obj <- Seurat::CreateSeuratObject(counts = pseudobulked_expr, assay = "scennep")
        SummarizedExperiment::assays(obj)[["scennep"]] <- pseudobulked_expr
    }
    return(obj)
}

extract_expression <- function(obj, assay, normalize_data) {
    # Map of class to function calls
    if (inherits(obj, "Seurat")) {
        layer <- ifelse(normalize_data, "data", "counts")
        expr <- Seurat::GetAssayData(obj, layer = layer)
    } else if (inherits(obj, "SingleCellExperiment")) {
        expr <- SummarizedExperiment::assay(obj, assay = assay)
    } else {
        stop("Unsupported object class")
    }

    return(expr)
}


build_snn <- function(obj, nn_count, npcs) {
    message("Building SNN graph with k = ", nn_count)
    snn_graph <- switch(
        class(obj)[1],
        "Seurat" = build_snn_seurat(obj, nn_count, npcs),
        "SingleCellExperiment" = build_snn_sce(obj, nn_count, npcs),
        stop("Unsupported object class")
    )
    return(snn_graph)
}

build_snn_seurat <- function(seu_obj, nn_count, npcs) {

    # detect_snn <- grepl("_snn", SeuratObject::Graphs(seu_obj))
    # if (!any(detect_snn)) {
    seu_obj <- Seurat::FindNeighbors(
        seu_obj,
        reduction = "pca",
        dims = 1:npcs,
        k.param = nn_count,
        compute.SNN = TRUE,
        prune.SNN = 0)
    slot <- SeuratObject::Graphs(seu_obj)[
        grepl("_snn", SeuratObject::Graphs(seu_obj))]
    # } else {
    #     message("Reusing existing SNN graph")
    #     slot <- SeuratObject::Graphs(seu_obj)[which(detect_snn)]
    # }

    # Extract SNN graph
    snn_graph <- SeuratObject::Graphs(seu_obj, slot = slot)
    snn_graph <- as(snn_graph, "dgCMatrix")
    return(snn_graph)
}


build_snn_sce <- function(sce, nn_count, npcs) {
    # Subset PCA
    SingleCellExperiment::reducedDim(sce, "PCA") <- SingleCellExperiment::reducedDim(sce, "PCA")[, 1:npcs, drop = FALSE]
    # Build the SNN graph
    snn_graph <- scran::buildSNNGraph(
        sce, 
        use.dimred = 'PCA', 
        k = nn_count,
        type = "jaccard")

    # Extract SNN graph
    snn_graph <- igraph::as_adjacency_matrix(snn_graph, attr = "weight")
    return(snn_graph)
}
