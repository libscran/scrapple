#' Normalize CRISPR counts
#'
#' Compute (log-)normalized expression values after performing scaling normalization of an CRISPR count matrix.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to CRISPR guides and columns correspond to cells.
#' @param size.factors Numeric vector of length equal to the number of columns of \code{x},
#' containing the size factor for each cell in \code{x}.
#' If \code{NULL}, this defaults to the column sums of the count matrix in \code{x}.
#' @param center Logical scalar indicating whether to center the \code{size.factors},
#' see \code{?\link[scrapper]{centerSizeFactors}} for more details.
#' @param block,mode Arguments passed to \code{\link[scrapper]{centerSizeFactors}}.
#' @param log,pseudo.count Arguments passed to \code{\link[scrapper]{normalizeCounts}}.
#' @param assay.type Integer or string specifying the assay of \code{x} with the count matrix.
#' @param output.name String containing the name of the assay to store the normalized matrix.
#' @param factor.name String containing the name of the \code{\link[SummarizedExperiment]{colData}} column in which to store the size factors in the output object.
#' If \code{NULL}, the size factors are not stored. 
#'
#' @return \code{x} is returned with a new assay containing the (log-)normalized matrix.
#' Size factors are also stored in the \code{\link[SummarizedExperiment]{colData}}.
#'
#' @author Aaron Lun
#'
#' @examples
#' sce <- altExp(getTestCrisprData.se("qc"), "CRISPR Guide Capture")
#' sce <- normalizeCrisprCounts.se(sce, size.factors=sce$sum)
#' assayNames(sce)
#' summary(sizeFactors(sce))
#'
#' @seealso
#' \code{\link[scrapper]{centerSizeFactors}} and \code{\link[scrapper]{normalizeCounts}}, from the \pkg{scrapper} package.
#'
#' @export
#' @importFrom SummarizedExperiment assay assay<-
#' @importFrom Matrix colSums
normalizeCrisprCounts.se <- function(
    x,
    size.factors = NULL,
    center = TRUE,
    block = NULL,
    mode = "lowest",
    log = TRUE,
    pseudo.count = 1,
    assay.type = "counts",
    output.name = "logcounts",
    factor.name = "sizeFactor"
) {
    y <- assay(x, assay.type)

    if (is.null(size.factors)) {
        size.factors <- colSums(y)
    }
    if (center) {
        size.factors <- scrapper::centerSizeFactors(size.factors, block=block, mode=mode)
    }

    assay(x, output.name) <- scrapper::normalizeCounts(y, size.factors=size.factors, log=log, pseudo.count=pseudo.count)
    if (!is.null(factor.name)) {
        colData(x)[[factor.name]] <- size.factors
    }

    x
}
