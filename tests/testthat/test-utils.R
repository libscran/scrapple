# library(testthat); library(scrapple); source("test-utils.R")

test_that(".call works as expected", {
    out <- scrapple:::.call(c, list(foo=1, bar=2), list(whee=3), list(stuff=4))
    expect_identical(out, c(foo=1, bar=2, whee=3, stuff=4))

    # Elements in more.args take precedence.
    out <- scrapple:::.call(c, list(foo=1, bar=2), list(whee=3), list(whee=2.5, stuff=4))
    expect_identical(out, c(foo=1, bar=2, whee=2.5, stuff=4))
})

test_that(".sanitize_altexp_assays works as expected", {
    out <- scrapple:::.sanitize_altexp_assays(1:2, c("foo", "bar", "whee"), "counts")
    expect_identical(out, c(foo="counts", bar="counts"))
    out <- scrapple:::.sanitize_altexp_assays(c("whee", "bar"), c("foo", "bar", "whee"), "counts")
    expect_identical(out, c(whee="counts", bar="counts"))

    # Handles duplicates.
    out <- scrapple:::.sanitize_altexp_assays(c(2, 1, 2), c("foo", "bar", "whee"), "counts")
    expect_identical(out, c(bar="counts", foo="counts"))

    # Works if a named vector is already supplied.
    out <- scrapple:::.sanitize_altexp_assays(c(foo="stuff", bar="whee"), c("foo", "bar", "whee"), "counts")
    expect_identical(out, c(foo="stuff", bar="whee"))

    # Still strips out duplicates, though.
    out <- scrapple:::.sanitize_altexp_assays(c(foo=2, bar=1, foo=3), c("foo", "bar", "whee"), "counts")
    expect_identical(out, c(foo=2, bar=1))
})

library(DelayedArray)
test_that(".add/get_transposed_reddim works as expected", {
    original <- SingleCellExperiment(list(counts=matrix(0, 10, 20)))
    rd <- matrix(runif(40), nrow=2)

    {
        sce <- original

        sce <- scrapple:::.add_transposed_reddim(sce, "foo", rd, delayed=FALSE)
        expect_identical(reducedDim(sce, "foo"), t(rd))
        expect_identical(scrapple:::.get_transposed_reddim(sce, "foo"), rd)

        sce <- scrapple:::.add_transposed_reddim(sce, "foo", rd, delayed=TRUE)
        expect_identical(reducedDim(sce, "foo"), t(DelayedArray(rd)))
        expect_identical(scrapple:::.get_transposed_reddim(sce, "foo"), rd)

        # Checking that transposition does indeed give us back a no-op on the original matrix.
        out <- t(reducedDim(sce, "foo"))
        expect_s4_class(out, "DelayedMatrix")
        expect_identical(out@seed, rd)

        # Works with other matrix types.
        reducedDim(sce, "stuff") <- t(Matrix(rd))
        expect_identical(scrapple:::.get_transposed_reddim(sce, "stuff"), rd)
    }

    # Trying again, with names.
    {
        sce <- original
        colnames(sce) <- LETTERS[1:20]
        named.rd <- rd
        colnames(named.rd) <- colnames(sce)

        sce <- scrapple:::.add_transposed_reddim(sce, "foo", rd, delayed=FALSE)
        expect_identical(scrapple:::.get_transposed_reddim(sce, "foo"), named.rd)

        sce <- scrapple:::.add_transposed_reddim(sce, "foo", rd, delayed=TRUE)
        expect_identical(scrapple:::.get_transposed_reddim(sce, "foo"), named.rd)

        reducedDim(sce, "stuff") <- t(Matrix(rd))
        expect_identical(scrapple:::.get_transposed_reddim(sce, "stuff"), named.rd)
    }

    # Get transposition works for altexps.
    {
        sce <- original
        rd2 <- matrix(runif(40), nrow=2)
        altExp(sce, "ADT") <- SingleCellExperiment(list(counts=matrix(0, 0, 20)), reducedDims=list(PCA=t(rd2)))
        expect_identical(scrapple:::.get_transposed_reddim(sce, c(ADT="PCA")), rd2)
    }
})
