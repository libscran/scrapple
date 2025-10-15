#' Quick quality control for CRISPR data
#'
#' Quickly compute quality control (QC) metrics, thresholds and filters from CRISPR data in a \link[SummarizedExperiment]{SummarizedExperiment}.
#' 
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to CRISPR guides and columns correspond to cells.
#' @param num.threads Arguments passed to \code{\link[scrapper]{computeCrisprQcMetrics}}.
#' @param num.mads,block Arguments passed to \code{\link[scrapper]{suggestCrisprQcThresholds}}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the CRISPR count matrix.
#' @param output.prefix String containing a prefix to add to the names of the \code{link[SummarizedExperiment]{colData}} columns containing the output statistics.
#' @param thresholds.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry containing the filtering thresholds.
#' If \code{NULL}, thresholds are not stored in the metadata.
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
    output.prefix = NULL,
    thresholds.name = "thresholds"
) {
    metrics <- scrapper::computeCrisprQcMetrics(assay(x, assay.type, withDimnames=FALSE), num.threads=num.threads)
    thresholds <- scrapper::suggestCrisprQcThresholds(metrics, block=block, num.mads=num.mads)
    keep <- scrapper::filterCrisprQcMetrics(thresholds, metrics, block=block)

    df <- formatComputeCrisprQcMetricsResult(metrics)
    df$keep <- keep
    colnames(df) <- paste0(output.prefix, colnames(df))
    colData(x) <- cbind(colData(x), df)

    if (!is.null(thresholds.name)) {
        metadata(x)[[thresholds.name]] <- thresholds
    }

    x
}

#' @export
#' @rdname quickCrisprQc.se
#' @importFrom S4Vectors DataFrame
formatComputeCrisprQcMetricsResult <- function(compute.res) {
    # Maybe we'll add something more later.
    DataFrame(compute.res)
}
