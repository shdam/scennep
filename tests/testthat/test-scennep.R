data("pbmc_small")
example_counts <- SeuratObject::GetAssayData(pbmc_small, layer = "counts")

test_that("scennep works with a matrix", {
    mat <- scennep(example_counts, nn_count = 5, return_S4 = FALSE)
    
    expect_type(mat, "double")
    expect_equal(dim(mat), dim(pbmc_small))
})

test_that("scennep works with Seurat", {
    seu <- scennep(pbmc_small, as = "seurat", nn_count = 5)
    
    expect_s4_class(seu, "Seurat")
    expect_true("scennep" %in% Seurat::Assays(seu), "Assay 'scennep' should be present")
    
})

test_that("scennep works with Bioc", {
    sce <- scennep(example_counts, as = "bioc", nn_count = 5)
    
    expect_s4_class(sce, "SingleCellExperiment")
    expect_true("scennep" %in% SummarizedExperiment::assayNames(sce), "Assay 'scennep' should be present")
})