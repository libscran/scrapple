#' Quick quality control for ADT data
#'
#' Quickly compute quality control (QC) metrics, thresholds and filters from ADT data in a \link[SummarizedExperiment]{SummarizedExperiment}.
#' 
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to antibody-derived tags (ADTs) and columns correspond to cells.
#' @param subsets,num.threads Arguments passed to \code{\link[scrapper]{computeAdtQcMetrics}}.
#' @param num.mads,block Arguments passed to \code{\link[scrapper]{suggestAdtQcThresholds}}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the ADT count matrix.
#' @param prefix String containing a prefix to append to the name of each column corresponding to a QC metric in the \code{link[SummarizedExperiment]{colData}}.
#' @param flatten Logical scalar indicating whether to flatten the subset proportions into separate columns of the \code{link[SummarizedExperiment]{colData}}.
#' If \code{FALSE}, the subset proportions are stored in a nested \link[S4Vectors]{DataFrame}.
#' @param compute.res List returned by \code{\link[scrapper]{computeAdtQcMetrics}}.
#' 
#' @return
#' For \code{quickAdtQc.se}, \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link[scrapper]{computeAdtQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
#'
#' For \code{formatComputeAdtQcMetricsResult}, a \link[S4Vectors]{DataFrame} is returned with the per-cell QC metrics.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link[scrapper]{computeAdtQcMetrics}},
#' \code{\link[scrapper]{suggestAdtQcThresholds}},
#' and \code{\link[scrapper]{filterAdtQcMetrics}}
#' from the \pkg{scrapper} package.
#'
#' @examples
#' sce <- altExp(getTestAdtData.se(), "ADT")
#' sce <- quickAdtQc.se(sce, subsets=list(igg=grepl("IgG", rownames(sce))))
#' colData(sce)[,c("sum", "detected", "igg.sum")]
#' metadata(sce)$thresholds
#' summary(sce$keep)
#' 
#' @export
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment assay colData colData<-
quickAdtQc.se <- function( 
    x,
    subsets,
    num.threads = 1,
    num.mads = 3,
    block = NULL,
    assay.type = "counts",
    prefix = NULL, 
    flatten = TRUE
) {
    metrics <- scrapper::computeAdtQcMetrics(assay(x, assay.type, withDimnames=FALSE), subsets, num.threads=num.threads)
    thresholds <- scrapper::suggestAdtQcThresholds(metrics, block=block, num.mads=num.mads)
    keep <- scrapper::filterAdtQcMetrics(thresholds, metrics, block=block)

    colData(x) <- cbind(colData(x), formatComputeAdtQcMetricsResult(metrics, prefix=prefix, flatten=flatten))
    colData(x)[[paste0(prefix, "keep")]] <- keep
    metadata(x)[[paste0(prefix, "thresholds")]] <- thresholds
    x
}

#' @export
#' @rdname quickAdtQc.se
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
formatComputeAdtQcMetricsResult <- function(compute.res, prefix = NULL, flatten = TRUE) {
    df <- DataFrame(sum=compute.res$sum, detected=compute.res$detected)

    if (flatten) {
        for (sub in names(compute.res$subsets)) {
            df[[paste0(sub, ".sum")]] <- compute.res$subsets[[sub]]
        }
    } else {
        tmp <- S4Vectors::make_zero_col_DFrame(nrow=nrow(cd))
        for (sub in names(compute.res$subsets)) {
            tmp[[paste0(sub)]] <- compute.res$subsets[[sub]]
        }
        df[["sum"]] <- tmp
    }

    colnames(df) <- paste0(prefix, colnames(df))
    df
}
