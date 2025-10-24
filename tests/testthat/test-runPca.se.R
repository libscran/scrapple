# library(testthat); library(scrapple); source("test-runPca.se.R")

set.seed(8888)
sce <- SingleCellExperiment(list(logcounts = matrix(rnorm(100000), 500, 200)))

test_that("runPca.se works as expected", {
    out <- runPca.se(sce, features=1:50) 
    expect_true("PCA" %in% reducedDimNames(out))
    expect_identical(nrow(metadata(out)$PCA$rotation), 50L)
    expect_identical(length(metadata(out)$PCA$variance.explained), 25L)

    # Works with features set to NULL.
    null <- runPca.se(sce[1:50,], features=NULL) 
    expect_identical(reducedDim(null), reducedDim(out))
    expect_identical(metadata(null)$PCA$rotation, metadata(out)$PCA$rotation)

    no.meta <- runPca.se(sce, features=1:10, meta.name=NULL)
    expect_null(metadata(no.meta)$PCA)
})
