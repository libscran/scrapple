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
#' @param raw Logical scalar indicating whether to return the QC results directly.
#' 
#' @return
#' For \code{quickAdtQc.se}:
#' \itemize{
#' \item If \code{raw=FALSE}, \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link[scrapper]{computeAdtQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
#' \item If \code{raw=TRUE}, a list is returned containing \code{metrics}, the result of \code{\link[scrapper]{computeAdtQcMetrics}};
#' \code{thresholds}, the result of \code{\link[scrapper]{suggestAdtQcThresholds}},
#' and \code{keep}, the result of \code{\link[scrapper]{FilterAdtQcMetrics}}.
#' }
#'
#' For \code{attachAdtQcMetrics.se}. \code{x} is returned with additional columns added to its \code{colData}.
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
#' # We can also manually execute many of the steps:
#' sce <- altExp(getTestAdtData.se(), "ADT")
#' res <- quickAdtQc.se(sce, subsets=list(igg=grepl("IgG", rownames(sce))), raw=TRUE)
#' sce <- attachAdtQcMetrics.se(sce, res$metrics, prefix="adt.qc.")
#' res$thresholds$detected <- 50 # manual override of the sum threshold
#' res$thresholds$subsets["igg"] <- 10 # same for the IgG sums
#' sce$adt.qc.keep <- scrapper::filterAdtQcMetrics(res$thresholds, res$metrics)
#' 
#' @export
#' @importFrom SummarizedExperiment assay colData<-
#' @importFrom SingleCellExperiment altExp
#' @importFrom S4Vectors metadata metadata<-
quickAdtQc.se <- function( 
    x,
    subsets = list(),
    num.threads = 1,
    num.mads = 3,
    block = NULL,
    assay.type = "counts",
    prefix = NULL, 
    flatten = TRUE,
    raw = FALSE
) {
    metrics <- scrapper::computeAdtQcMetrics(assay(x, assay.type, withDimnames=FALSE), subsets, num.threads=num.threads)
    thresholds <- scrapper::suggestAdtQcThresholds(metrics, block=block, num.mads=num.mads)
    keep <- scrapper::filterAdtQcMetrics(thresholds, metrics, block=block)

    if (raw) {
        return(list(metrics=metrics, thresholds=thresholds, keep=keep))
    }

    x <- attachAdtQcMetrics.se(x, metrics, prefix=prefix, flatten=flatten)
    colData(x)[[paste0(prefix, "keep")]] <- keep
    metadata(x)[[paste0(prefix, "thresholds")]] <- thresholds
    x
}

#' @export
#' @rdname quickAdtQc.se
#' @importFrom SummarizedExperiment colData colData<-
attachAdtQcMetrics.se <- function(x, metrics, prefix = "", flatten = TRUE) {
    cd <- colData(x)
    cd[[paste0(prefix, "sum")]] <- metrics$sum
    cd[[paste0(prefix, "detected")]] <- metrics$detected

    if (flatten) {
        for (sub in names(metrics$subsets)) {
            cd[[paste0(prefix, sub, ".sum")]] <- metrics$subsets[[sub]]
        }
    } else {
        tmp <- S4Vectors::make_zero_col_DFrame(nrow=nrow(cd))
        for (sub in names(metrics$subsets)) {
            tmp[[paste0(sub)]] <- metrics$subsets[[sub]]
        }
        cd[[paste0(prefix, "sum")]] <- tmp
    }

    colData(x) <- cd
    x
}
