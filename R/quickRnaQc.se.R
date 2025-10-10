#' Quick quality control for RNA data
#'
#' Quickly compute quality control (QC) metrics, thresholds and filters from RNA data in a \link[SummarizedExperiment]{SummarizedExperiment}.
#' 
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @inheritParams scrapper::computeRnaQcMetrics
#' @inheritParams scrapper::suggestRnaQcThresholds
#' @param assay.type Integer or string specifying the assay of \code{x} containing the RNA count matrix.
#' @param altexp.proportions Character vector containing the names of alternative experiments for which to compute proportions relative to the RNA total.
#' For each alternative experiment, the proportion is defined \deqn{X/(X+Y)} where \deqn{X} is the alternative experiment's total and \deqn{Y} is the RNA total.
#' The count matrix in each alternative experiment should be stored in the assay specified by \code{assay.type}.
#' This parameter is typically used to refer to alternative experiments holding spike-in data. 
#' @param prefix String containing a prefix to append to the name of each column corresponding to a QC metric in the \code{link[SummarizedExperiment]{colData}}.
#' @param flatten Logical scalar indicating whether to flatten the subset proportions into separate columns of the \code{link[SummarizedExperiment]{colData}}.
#' If \code{FALSE}, the subset proportions are stored in a nested \link[S4Vectors]{DataFrame}.
#' @param raw Logical scalar indicating whether to return the QC results directly.
#' 
#' @return
#' For \code{quickRnaQc.se}:
#' \itemize{
#' \item If \code{raw=FALSE}, \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link[scrapper]{computeRnaQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
#' \item If \code{raw=TRUE}, a list is returned containing \code{metrics}, the result of \code{\link[scrapper]{computeRnaQcMetrics}};
#' \code{thresholds}, the result of \code{\link[scrapper]{suggestRnaQcThresholds}},
#' and \code{keep}, the result of \code{\link[scrapper]{FilterRnaQcMetrics}}.
#' }
#'
#' For \code{attachRnaQcMetrics.se}. \code{x} is returned with additional columns added to its \code{colData}.
#'
#' @author Aaron Lun
#'
#' @seealso
#' \code{\link[scrapper]{computeRnaQcMetrics}},
#' \code{\link[scrapper]{suggestRnaQcThresholds}},
#' and \code{\link[scrapper]{filterRnaQcMetrics}}
#' from the \pkg{scrapper} package.
#'
#' @examples
#' sce <- getTestRnaData.se()
#' sce <- quickRnaQc.se(sce, subsets=list(mito=grepl("^mt", rownames(sce))))
#' colData(sce)[,c("sum", "detected", "mito.proportion")]
#' metadata(sce)$thresholds
#' summary(sce$keep)
#'
#' # Computing spike-in proportions, if available.
#' sce <- getTestRnaData.se()
#' sce <- quickRnaQc.se(
#'    sce,
#'    subsets=list(mito=grepl("^mt", rownames(sce))),
#'    altexp.proportions="ERCC"
#' )
#' colData(sce)[,c("sum", "detected", "mito.proportion", "ERCC.proportion")]
#'
#' # We can also manually execute many of the steps:
#' sce <- getTestRnaData.se()
#' res <- quickRnaQc.se(sce, subsets=list(mito=grepl("^mt", rownames(sce))), raw=TRUE)
#' sce <- attachRnaQcMetrics.se(sce, res$metrics, prefix="rna.qc.")
#' res$thresholds$sum <- 1000 # manual override of the sum threshold
#' res$thresholds$subsets["mito"] <- 0.05 # same for the mitochondrial proportions
#' sce$rna.qc.keep <- scrapper::filterRnaQcMetrics(res$thresholds, res$metrics)
#' 
#' @export
#' @importFrom SummarizedExperiment assay colData<-
#' @importFrom SingleCellExperiment altExp
#' @importFrom S4Vectors metadata metadata<-
quickRnaQc.se <- function( 
    x,
    subsets = list(),
    block = NULL,
    num.mads = 3,
    num.threads = 1,
    assay.type = "counts",
    altexp.proportions = NULL,
    prefix = NULL, 
    flatten = TRUE,
    raw = FALSE
) {
    metrics <- scrapper::computeRnaQcMetrics(assay(x, assay.type, withDimnames=FALSE), subsets, num.threads=num.threads)

    # Adding more proportions from the alternative experiments, mostly for spike-ins.
    if (!is.null(altexp.proportions)) {
        altexp.names <- names(altexp.proportions)
        if (!is.null(altexp.names)) {
            altexp.names <- altexp.proportions
        }
        for (i in seq_along(altexp.proportions)) {
            alt.assay <- assay(altExp(x, altexp.proportions[i]), withDimnames=FALSE)
            alt.metrics <- scrapper::computeRnaQcMetrics(alt.assay, subsets=list(), num.threads=num.threads)
            metrics$subsets[[altexp.names[i]]] <- alt.metrics$sum/(metrics$sum + alt.metrics$sum)
        }
    }

    thresholds <- scrapper::suggestRnaQcThresholds(metrics, block=block, num.mads=num.mads)
    keep <- scrapper::filterRnaQcMetrics(thresholds, metrics, block=block)

    if (raw) {
        return(list(metrics=metrics, thresholds=thresholds, keep=keep))
    }

    x <- attachRnaQcMetrics.se(x, metrics, prefix=prefix, flatten=flatten)
    colData(x)[[paste0(prefix, "keep")]] <- keep
    metadata(x)[[paste0(prefix, "thresholds")]] <- thresholds
    x
}

#' @export
#' @rdname quickRnaQc.se
#' @importFrom SummarizedExperiment colData colData<-
attachRnaQcMetrics.se <- function(x, metrics, prefix = "", flatten = TRUE) {
    cd <- colData(x)
    cd[[paste0(prefix, "sum")]] <- metrics$sum
    cd[[paste0(prefix, "detected")]] <- metrics$detected

    if (flatten) {
        for (sub in names(metrics$subsets)) {
            cd[[paste0(prefix, sub, ".proportion")]] <- metrics$subsets[[sub]]
        }
    } else {
        tmp <- S4Vectors::make_zero_col_DFrame(nrow=nrow(cd))
        for (sub in names(metrics$subsets)) {
            tmp[[paste0(sub)]] <- metrics$subsets[[sub]]
        }
        cd[[paste0(prefix, "proportion")]] <- tmp
    }

    colData(x) <- cd
    x
}
