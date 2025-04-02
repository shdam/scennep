
test_that("scennep works with Seurat objects", {
  # Create a small example dataset
  set.seed(42)
  data <- matrix(rpois(1000, lambda = 5), ncol = 100, nrow = 10)
  rownames(data) <- paste0("Gene", 1:10)
  colnames(data) <- paste0("Cell", 1:100)
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = data)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  
  # Run scennep
  result <- scennep(seurat_obj, nn_count = 10, return_seurat = TRUE)
  
  # Tests
  expect_s4_class(result, "Seurat")
  expect_true("scennep" %in% Assays(result))
  expect_equal(dim(GetAssayData(result, assay = "scennep")), dim(data))
  expect_true(all(colnames(GetAssayData(result, assay = "scennep")) == colnames(data)))
})

test_that("scennep works with SingleCellExperiment objects", {
  # Create a small example dataset
  set.seed(42)
  data <- matrix(rpois(1000, lambda = 5), ncol = 100, nrow = 10)
  rownames(data) <- paste0("Gene", 1:10)
  colnames(data) <- paste0("Cell", 1:100)
  
  # Create a SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = data))
  sce <- scuttle::logNormCounts(sce)
  
  # Run scennep
  result <- scennep(sce, nn_count = 10, as = "bioc", assay = "counts")
  
  # Tests
  expect_s4_class(result, "SingleCellExperiment")
  expect_true("scennep" %in% assayNames(result))
  expect_equal(dim(assay(result, "scennep")), dim(data))
  expect_true(all(colnames(assay(result, "scennep")) == colnames(data)))
})

test_that("scennep handles edge cases", {
  # Create a small example dataset
  set.seed(42)
  data <- matrix(rpois(100, lambda = 5), ncol = 10, nrow = 10)
  rownames(data) <- paste0("Gene", 1:10)
  colnames(data) <- paste0("Cell", 1:10)
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = data)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  
  # Test with nn_count larger than number of cells
  result <- scennep(seurat_obj, nn_count = 20, return_seurat = TRUE)
  expect_s4_class(result, "Seurat")
  expect_true("scennep" %in% Assays(result))
  
  # Test with single-cell input
  single_cell_data <- matrix(rpois(10, lambda = 5), ncol = 1, nrow = 10)
  single_cell_obj <- CreateSeuratObject(counts = single_cell_data)
  expect_error(scennep(single_cell_obj, nn_count = 10), "Not enough cells for pseudobulking")
})

test_that("scennep respects markers parameter", {
  # Create a small example dataset
  set.seed(42)
  data <- matrix(rpois(1000, lambda = 5), ncol = 100, nrow = 10)
  rownames(data) <- paste0("Gene", 1:10)
  colnames(data) <- paste0("Cell", 1:100)
  
  # Create a Seurat object
  seurat_obj <- CreateSeuratObject(counts = data)
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(seurat_obj)
  seurat_obj <- ScaleData(seurat_obj)
  seurat_obj <- RunPCA(seurat_obj)
  
  # Run scennep with specific markers
  markers <- c("Gene1", "Gene3", "Gene5")
  result <- scennep(seurat_obj, nn_count = 10, markers = markers, return_seurat = TRUE)
  
  # Tests
  expect_equal(rownames(GetAssayData(result, assay = "scennep")), markers)
})
