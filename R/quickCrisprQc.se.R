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
#' @param compute.res List returned by \code{\link[scrapper]{computeCrisprQcMetrics}}.
#' 
#' @return
#' For \code{quickCrisprQc.se}, \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link[scrapper]{computeCrisprQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
#'
#' For \code{formatComputeCrisprQcMetricsResult}, a \link[S4Vectors]{DataFrame} is returned with the per-cell QC metrics.
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
#' @export
#' @importFrom SummarizedExperiment assay colData colData<-
#' @importFrom S4Vectors metadata metadata<-
quickCrisprQc.se <- function( 
    x,
    num.threads = 1,
    num.mads = 3,
    block = NULL,
    assay.type = "counts",
    prefix = NULL, 
    flatten = TRUE
) {
    metrics <- scrapper::computeCrisprQcMetrics(assay(x, assay.type, withDimnames=FALSE), num.threads=num.threads)
    thresholds <- scrapper::suggestCrisprQcThresholds(metrics, block=block, num.mads=num.mads)
    keep <- scrapper::filterCrisprQcMetrics(thresholds, metrics, block=block)

    colData(x) <- cbind(colData(x), formatComputeCrisprQcMetricsResult(metrics, prefix=prefix, flatten=flatten))
    colData(x)[[paste0(prefix, "keep")]] <- keep
    metadata(x)[[paste0(prefix, "thresholds")]] <- thresholds
    x
}

#' @export
#' @rdname quickCrisprQc.se
#' @importFrom S4Vectors DataFrame
formatComputeCrisprQcMetricsResult <- function(compute.res, prefix = NULL, flatten = TRUE) {
    df <- DataFrame(compute.res)
    colnames(df) <- paste0(prefix, colnames(df))
    df
}
