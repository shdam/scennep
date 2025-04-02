#' Pseudobulk Neighbors for Single-Cell Expression Data
#'
#' This function takes a Seurat object with a precomputed SNN graph and computes pseudobulk
#' expression profiles for each cell based on its neighbors. The pseudobulk profiles are
#' created by summing the expression values of a cell and its top nearest neighbors, specified
#' by the `nn_count`. If a cell has fewer neighbors than `nn_count`, all available neighbors
#' are used. The result is added as a new assay in the Seurat object.
#'
#' @param seurat_obj A Seurat object containing single-cell expression data with an
#'   existing SNN graph stored in the "SCT_snn" slot.
#' @param nn_count The number of nearest neighbors to consider for pseudobulking. If a
#'   cell has fewer than this number of neighbors, all are used. Default is 100.
#' @param mc.cores The number of cores to use for parallel processing. By default, it uses
#'   one less than the total number of cores available to prevent the system from locking up.
#'
#' @return The input Seurat object with an additional assay named "SCT_bulk", which contains
#'   the pseudobulked expression data.
#'
#' @examples
#' # Assuming `pbmc` is a Seurat object with SCTransform normalization and an SNN graph computed
#' pbmc <- pseudobulk_neighbors(pbmc, nn_count = 100, mc.cores = 8)
#'
#' @export
pseudobulk_neighbors <- function(
    seurat_obj,
    nn_count = 100,
    mc.cores = parallel::detectCores() - 1) {
  if (!"SCT_snn" %in% Graphs(seurat_obj)) {
    message("Computing SNN graph")
    seurat_obj <- FindNeighbors(seurat_obj, reduction = "pca", compute.SNN = TRUE)
  }

  snn_graph <- Graphs(seurat_obj, slot = "SCT_snn")
  if (!inherits(snn_graph, "dgCMatrix")) {
    snn_graph <- as(snn_graph, "dgCMatrix")
  }

  num_cells <- ncol(snn_graph)
  expr <- GetAssayData(seurat_obj, assay = "SCT", layer = "data")

  message("Pseudo-bulking each cell with 100 neighbors")
  pseudobulked_expr <- pbmcapply::pbmclapply(1:num_cells, function(cell_id) {
    neighbors <- which(snn_graph[, cell_id] > 0)
    if (length(neighbors) > nn_count) {
      neighbors <- sample(neighbors, nn_count)
    }
    pseudobulked <- rowSums(expr[, c(cell_id, neighbors)])
    return(pseudobulked)
  }, mc.cores = mc.cores)

  message("Adding pseudo-bulked expression data to 'SCT_bulk'")
  rm(expr)
  pseudobulked_expr <- do.call(cbind, pseudobulked_expr[[1]])
  colnames(pseudobulked_expr) <- colnames(seurat_obj)

  pseudobulked_expr <- CreateAssayObject(counts = pseudobulked_expr, key = "sctbulk_")
  seurat_obj[["SCT_bulk"]] <- pseudobulked_expr

  return(seurat_obj)
}
