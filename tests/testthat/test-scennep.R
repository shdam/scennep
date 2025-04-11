data("example_counts")


test_that("scennep works with Seurat", {
    seu <- scennep(example_counts, as = "seu")
    
    expect_s4_class(seu, "Seurat")
    expect_true("scennep" %in% Seurat::Assays(seu), "Assay 'scennep' should be present")
    
})

test_that("scennep works with Bioc", {
    sce <- scennep(example_counts, as = "bioc")
    
    expect_s4_class(sce, "SingleCellExperiment")
    expect_true("scennep" %in% SummarizedExperiment::assayNames(sce), "Assay 'scennep' should be present")
})