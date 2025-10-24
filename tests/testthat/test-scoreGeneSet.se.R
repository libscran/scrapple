# library(testthat); library(scrapple); source("test-scoreGeneSet.se.R")

set.seed(1111)
sce <- SummarizedExperiment(list(logcounts = matrix(runif(2000), 200, 10)))

test_that("scoreGeneSet.se works as expected", {
    out <- scoreGeneSet.se(sce, set=1:10)
    expect_identical(length(out$scores), ncol(sce))
    expect_identical(length(out$weights), 10L)

    # Respects names in the output object.
    rownames(sce) <- paste0("GENE_", seq_len(nrow(sce)))
    out2 <- scoreGeneSet.se(sce, set=10:1)
    expect_equal(out$scores, out2$scores)
    expect_identical(names(out2$weights), rownames(sce)[1:10])
    expect_identical(unname(out2$weights), out$weights)
})
