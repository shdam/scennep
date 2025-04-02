#' Pseudobulk Neighbors for Single-Cell Expression Data
#'
#' This function takes a Seurat or SCE object and computes pseudobulk
#' expression profiles for each cell based on its neighbors. The pseudobulk profiles are
#' created by summing the expression values of a cell and its top nearest neighbors, specified
#' by the `nn_count`. If a cell has fewer neighbors than `nn_count`, all available neighbors
#' are used. The result is added as a new assay in the Seurat object.
#'
#' @param obj A Seurat or SCE object or matrix containing single-cell expression data.
#' @param nn_count The number of nearest neighbors to consider for pseudobulking. If a
#'   cell has fewer than this number of neighbors, all are used. Default is 100.
#' @param nn_cutoff The strength cutoff for shared neighbors.
#' @param FUN Aggregation function. Defaults to rowSums.
#' @param npcs Number of PCs to use for the neighborhood graph.
#' @param markers Subset to markers of interest.
#' @param as If you provide a matrix, decide if you want to use Seurat or Bioconductor methods.
#' @param assay The name of the assay in the SingleCellExperiment object.
#' @param mc.cores The number of cores to use for parallel processing. By default, it uses
#'   one less than the total number of cores available to prevent the system from locking up.
#' @param pb Logical. If TRUE, a progress bar is printed.
#' @param return_seurat Logical. If TRUE, the Seurat object is returned with an nnbulk assay, but this is at the cost of memmory efficiency.
#' I recommend setting it to FALSE and adding the assay manually.
#'
#' @return The input Seurat object with an additional assay named "SCT_bulk", which contains
#'   the pseudobulked expression data.
#'
#' @examples
#' # Assuming `pbmc` is a Seurat object with SCTransform normalization and an SNN graph computed
#' pbmc <- pseudobulk_neighbors(pbmc, nn_count = 100, mc.cores = 8, return_seurat = TRUE)
#'
#' @export
pseudobulk_neighbors <- function(
    obj,
    nn_count = 100,
    nn_cutoff = 1/2,
    npcs = 15,
    FUN = rowSums,
    markers = NULL,
    as = c("seurat", "bioc"),
    return_seurat = TRUE,
    normalize_data = TRUE,
    assay = "SCT",
    pb = TRUE,
    mc.cores = 1) {

  if(pb) {
    cyCombine:::missing_package(package = "pbmcapply")
    APPLY <- pbmcapply::pbmclapply
  } else {
    APPLY <- parallel::mclapply
  }
  formals(APPLY)$mc.cores <- mc.cores

  stopifnot("Valid object types are 'Seurat', 'matrix', 'SingleCellExperiment'"
            = class(obj)[1] %in%
              c("Seurat", "matrix", "dgCMatrix", "SingleCellExperiment"))
  as <- match.arg(as)
  if (class(obj)[1] != "Seurat") {
    return_seurat <- FALSE
  }
  if (class(obj)[1] %in% c("matrix", "dgCMatrix")) {
    if (as == "seurat") {
      obj <- Seurat::CreateSeuratObject(counts = obj)
      obj <- FindVariableFeatures(obj)
      # if (!is.null(markers)) Seurat::VariableFeatures(obj) <- markers
    } else {
      obj <- SingleCellExperiment::SingleCellExperiment(obj)
      SummarizedExperiment::assayNames(obj) <- assay
    }
  }

  # Build PCA
  if (class(obj)[1] == "Seurat") {
    if (!"SCT" %in% Seurat::Assays(obj)) {
      # if (!"data" %in% Layers(obj)) obj <- NormalizeData(obj)
      # if (!"scale.data" %in% Layers(obj)) obj <- ScaleData(obj)

      obj <- Seurat::SCTransform(obj)
    } else {
      Seurat::DefaultAssay(obj) <- "SCT"
      if ("nnbulk" %in% Seurat::Assays(obj)){
        warning("Assay 'nnbulk' already exists! It will be overwritten.\n")
        obj[["nnbulk"]] <- NULL
      }
    }
    # DefaultAssay(obj) <- "SCT"
    if (length(SeuratObject::VariableFeatures(obj)) == 0) obj <- Seurat::FindVariableFeatures(obj)
    if(!"scale.data" %in% SeuratObject::Layers(obj)) obj <-  Seurat::ScaleData(obj)
    if (!"pca" %in% Seurat::Reductions(obj)) obj <- Seurat::RunPCA(obj, npcs = npcs)
  } else {
    if (assay == "counts" & normalize_data) {
      obj <- scuttle::logNormCounts(obj)
      assay <- "logcounts"
    }
    obj <- scater::runPCA(obj, ncomponents = npcs, scale = TRUE, exprs_values = assay)
  }


  # Build SNN graph
  snn_graph <- build_snn(obj, nn_count, npcs)


  # Extract expression data
  expr <- extract_expression(obj, assay, normalize_data)

  if (is.null(markers)) markers <- rownames(exprs)
  expr <- expr[markers, ]
  cols <- colnames(obj)

  # More memmory efficient to rm Seurat object
  if (!return_seurat) {
    rm(obj)
    gc()
  }
  num_cells <- ncol(snn_graph)
  message("Pseudo-bulking each cell with its ", nn_count, " neighbors")
  pseudobulked_expr <- APPLY(seq_len(num_cells), function(cell_id) {
    neighbors <- which(snn_graph[, cell_id] > nn_cutoff)

    if (length(neighbors) > nn_count) {
      # Order neighbors by SNN strength
      neighbor_weights <- snn_graph[neighbors, cell_id]
      top_neighbors <- order(neighbor_weights, decreasing = TRUE)[seq_len(nn_count)]
      neighbors <- neighbors[top_neighbors]
    }
    if (length(neighbors) == 1) neighbors <- c(neighbors, neighbors)
    pseudobulked <- FUN(expr[markers, neighbors])
    return(pseudobulked)
  })

  rm(expr)
  pseudobulked_expr <- do.call(cbind, pseudobulked_expr)
  colnames(pseudobulked_expr) <- cols

  if (!return_seurat) return(pseudobulked_expr)

  message("Adding pseudo-bulked expression data to assay 'nnbulk'")
  pseudobulked_expr <- Seurat::CreateAssayObject(data = pseudobulked_expr, key = "nnbulk_")
  obj[["nnbulk"]] <- pseudobulked_expr
  Seurat::DefaultAssay(obj) <- "nnbulk"

  return(obj)
}

extract_expression <- function(obj, assay, normalize_data) {
  # Map of class to function calls
  if (class(obj)[1] == "Seurat") {
    # assay <- ifelse(normalize_data, "SCT", "RNA")
    layer <- ifelse(normalize_data, "data", "counts")
    expr <- Seurat::GetAssayData(obj, assay = assay, layer = layer)
  } else if (class(obj)[1] == "SingleCellExperiment") {
    expr <- SummarizedExperiment::assay(obj, assay = assay)
  } else {
    stop("Unsupported object class")
  }

  # expr <- switch(class(obj)[1],
  #                "Seurat" = Seurat::GetAssayData(obj, assay = "SCT", layer = "data"),
  #                "SingleCellExperiment" = SummarizedExperiment::assay(obj, assay = assay),
  #                stop("Unsupported object class")
  # )
  return(expr)
}


build_snn <- function(obj, nn_count, npcs) {
  snn_graph <- switch(class(obj)[1],
                      "Seurat" = build_snn_seurat(obj, nn_count, npcs),
                      "SingleCellExperiment" = build_snn_sce(obj, nn_count, npcs),
                      stop("Unsupported object class")
  )
  return(snn_graph)
}

build_snn_seurat <- function(seu_obj, nn_count, npcs) {

  if (!"SCT_snn" %in% SeuratObject::Graphs(seu_obj)) {
    message("Building SNN graph")
    seu_obj <- Seurat::FindNeighbors(
      seu_obj,
      reduction = "pca",
      # dims = npcs,
      # k.param = nn_count,
      compute.SNN = TRUE)
  }

  # Extract SNN graph
  snn_graph <- SeuratObject::Graphs(seu_obj, slot = "SCT_snn")
  snn_graph <- as(snn_graph, "dgCMatrix")
  return(snn_graph)
}


build_snn_sce <- function(sce, nn_count, assay, npcs) {
  # Build the SNN graph
  snn_graph <- scran::buildSNNGraph(sce, use.dimred = 'PCA', k = nn_count, type = "jaccard")

  # Extract the SNN graph
  snn_graph <- igraph::as_adj(snn_graph)
  return(snn_graph)
}
