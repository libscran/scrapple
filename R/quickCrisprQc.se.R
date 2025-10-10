#' Quick quality control for CRISPR data
#'
#' Quickly compute quality control (QC) metrics, thresholds and filters from CRISPR data in a \link[SummarizedExperiment]{SummarizedExperiment}.
#' 
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to CRISPR guides and columns correspond to cells.
#' @param num.threads Arguments passed to \code{\link[scrapper]{computeCrisprQcMetrics}}.
#' @param num.mads,block Arguments passed to \code{\link[scrapper]{suggestCrisprQcThresholds}}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the CRISPR count matrix.
#' @param prefix String containing a prefix to append to the name of each column corresponding to a QC metric in the \code{link[SummarizedExperiment]{colData}}.
#' @param flatten Logical scalar indicating whether to flatten the subset proportions into separate columns of the \code{link[SummarizedExperiment]{colData}}.
#' If \code{FALSE}, the subset proportions are stored in a nested \link[S4Vectors]{DataFrame}.
#' @param raw Logical scalar indicating whether to return the QC results directly.
#' 
#' @return
#' For \code{quickCrisprQc.se}:
#' \itemize{
#' \item If \code{raw=FALSE}, \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link[scrapper]{computeCrisprQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
#' \item If \code{raw=TRUE}, a list is returned containing \code{metrics}, the result of \code{\link[scrapper]{computeCrisprQcMetrics}};
#' \code{thresholds}, the result of \code{\link[scrapper]{suggestCrisprQcThresholds}},
#' and \code{keep}, the result of \code{\link[scrapper]{FilterCrisprQcMetrics}}.
#' }
#'
#' For \code{attachCrisprQcMetrics.se}. \code{x} is returned with additional columns added to its \code{colData}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link[scrapper]{computeCrisprQcMetrics}},
#' \code{\link[scrapper]{suggestCrisprQcThresholds}},
#' and \code{\link[scrapper]{filterCrisprQcMetrics}}
#' from the \pkg{scrapper} package.
#'
#' @examples
#' sce <- altExp(getTestCrisprData.se(), "CRISPR Guide Capture")
#' sce <- quickCrisprQc.se(sce)
#' colData(sce)[,c("sum", "detected", "max.value", "max.index")]
#' metadata(sce)$thresholds
#' summary(sce$keep)
#'
#' # We can also manually execute many of the steps:
#' sce <- altExp(getTestCrisprData.se(), "CRISPR Guide Capture")
#' res <- quickCrisprQc.se(sce, raw=TRUE)
#' sce <- attachCrisprQcMetrics.se(sce, res$metrics, prefix="crispr.qc.")
#' res$thresholds$max.value <- 50 # manual override of the max value threshold
#' sce$crispr.qc.keep <- scrapper::filterCrisprQcMetrics(res$thresholds, res$metrics)
#' 
#' @export
#' @importFrom SummarizedExperiment assay colData<-
#' @importFrom SingleCellExperiment altExp
#' @importFrom S4Vectors metadata metadata<-
quickCrisprQc.se <- function( 
    x,
    num.threads = 1,
    num.mads = 3,
    block = NULL,
    assay.type = "counts",
    prefix = NULL, 
    flatten = TRUE,
    raw = FALSE
) {
    metrics <- scrapper::computeCrisprQcMetrics(assay(x, assay.type, withDimnames=FALSE), num.threads=num.threads)
    thresholds <- scrapper::suggestCrisprQcThresholds(metrics, block=block, num.mads=num.mads)
    keep <- scrapper::filterCrisprQcMetrics(thresholds, metrics, block=block)

    if (raw) {
        return(list(metrics=metrics, thresholds=thresholds, keep=keep))
    }

    x <- attachCrisprQcMetrics.se(x, metrics, prefix=prefix, flatten=flatten)
    colData(x)[[paste0(prefix, "keep")]] <- keep
    metadata(x)[[paste0(prefix, "thresholds")]] <- thresholds
    x
}

#' @export
#' @rdname quickCrisprQc.se
#' @importFrom SummarizedExperiment colData colData<-
attachCrisprQcMetrics.se <- function(x, metrics, prefix = "", flatten = TRUE) {
    cd <- colData(x)
    cd[[paste0(prefix, "sum")]] <- metrics$sum
    cd[[paste0(prefix, "detected")]] <- metrics$detected
    cd[[paste0(prefix, "max.value")]] <- metrics$max.value
    cd[[paste0(prefix, "max.index")]] <- metrics$max.index
    colData(x) <- cd
    x
}
