#' Normalize ADT counts
#'
#' Compute (log-)normalized expression values after performing scaling normalization of an ADT count matrix.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to antibody-derived tags (ADTs) and columns correspond to cells.
#' @param size.factors Numeric vector of length equal to the number of columns of \code{x},
#' containing the size factor for each cell in \code{x}.
#' If \code{NULL}, this defaults to the output of \code{\link[scrapper]{computeClrm1Factors}}.
#' @param num.threads Arguments passed to \code{\link[scrapper]{computeClrm1Factors}}.
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
#' sce <- altExp(getTestAdtData.se("qc"), "ADT")
#' sce <- normalizeAdtCounts.se(sce)
#' assayNames(sce)
#' summary(sizeFactors(sce))
#'
#' @seealso
#' \code{\link[scrapper]{computeClrm1Factors}},
#' \code{\link[scrapper]{centerSizeFactors}}
#' and \code{\link[scrapper]{normalizeCounts}}, from the \pkg{scrapper} package.
#'
#' @export
#' @importFrom SummarizedExperiment assay assay<- colData colData<-
normalizeAdtCounts.se <- function(
    x,
    size.factors = NULL,
    num.threads = 1,
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
        size.factors <- scrapper::computeClrm1Factors(y, num.threads=num.threads)
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
