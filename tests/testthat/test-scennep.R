data("example_counts")

test_that("scennep works with a matrix", {
    sen <- scennep(example_counts, nn_count = 5, return_S4 = FALSE)
    
    expect_type(sen, "double")
    expect_equal(dim(sen), dim(example_counts))
})

test_that("scennep works with Seurat", {
    seu <- scennep(example_counts, as = "seurat", nn_count = 5)
    
    expect_s4_class(seu, "Seurat")
    expect_true("scennep" %in% Seurat::Assays(seu), "Assay 'scennep' should be present")
    
})

test_that("scennep works with Bioc", {
    sce <- scennep(example_counts, as = "bioc", nn_count = 5)
    
    expect_s4_class(sce, "SingleCellExperiment")
    expect_true("scennep" %in% SummarizedExperiment::assayNames(sce), "Assay 'scennep' should be present")
})