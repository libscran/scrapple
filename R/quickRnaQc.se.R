#' Quick quality control for RNA data
#'
#' Quickly compute quality control (QC) metrics, thresholds and filters from RNA data in a \link[SummarizedExperiment]{SummarizedExperiment}.
#' 
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param subsets,num.threads Arguments passed to \code{\link[scrapper]{computeRnaQcMetrics}}.
#' @param num.mads,block Arguments passed to \code{\link[scrapper]{suggestRnaQcThresholds}}.
#' @param altexp.proportions Unnamed integer or character vector containing the indices/names of alternative experiments for which to compute QC metrics, see Details.
#' The assay to use from each alternative experiment is determined by \code{assay.type}.
#'
#' Alternatively, a named integer or character vector.
#' Each name specifies an alternative experiment while each value is the index/name of the assay to use from that experiment.
#' 
#' Only relevant if \code{x} is a \link[SingleCellExperiment]{SingleCellExperiment}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the RNA count matrix.
#' @param output.prefix String containing a prefix to add to the names of the \code{link[SummarizedExperiment]{colData}} columns containing the output statistics.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry containing the additional outputs such as the filtering thresholds.
#' If \code{NULL}, additional outputs are not reported. 
#' @param flatten Logical scalar indicating whether to flatten the subset proportions into separate columns of the \code{link[SummarizedExperiment]{colData}}.
#' If \code{FALSE}, the subset proportions are stored in a nested \link[S4Vectors]{DataFrame}.
#' @param compute.res List returned by \code{\link[scrapper]{computeRnaQcMetrics}}.
#' 
#' @details
#' \code{altexp.proportions} is typically used to refer to alternative experiments holding spike-in data.
#' For each alternative experiment, the proportion is defined \deqn{X/(X+Y)} where \deqn{X} is the alternative experiment's total and \deqn{Y} is the RNA total.
#' These proportions will be used for filtering in the same manner as the proportions computed from \code{subsets}.
#'
#' @return
#' For \code{quickRnaQc.se}, \code{x} is returned with additional columns added to its \code{\link[SummarizedExperiment]{colData}}.
#' Each column contains per-cell values for one of the QC metrics, see \code{\link[scrapper]{computeRnaQcMetrics}} for details.
#' The suggested thresholds are stored as a list in \code{\link[S4Vectors]{metadata}}.
#' The \code{colData} also contains a \code{keep} column, specifying which cells are to be retained.
#' If \code{altexp.proportions} is provided, QC metrics are added to the \code{colData} of the specified alternative experiments in the output object.
#'
#' For \code{computeRnaQcMetricsWithAltExps}, a list is returned containing:
#' \itemize{
#' \item \code{main}, the result of calling \code{\link[scrapper]{computeRnaQcMetrics}} on the RNA count matrix in \code{x}.
#' The proportion of counts in each alternative experiment is added to the \code{subsets}.
#' \item \code{altexp}, a named list of length equal to \code{altexp.proportions}.
#' Each inner list is the result of calling \code{computeRnaQcMetrics} on the RNA count matrix of the corresponding alternative experiment of \code{x}.
#' }
#'
#' For \code{formatComputeRnaQcMetricsResult}, a \link[S4Vectors]{DataFrame} is returned containing the per-cell QC metrics.
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
#' metadata(sce)$qc$thresholds
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
#' colData(altExp(sce, "ERCC"))[,c("sum", "detected")]
#'
#' @export
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment altExp altExp<- altExpNames
quickRnaQc.se <- function( 
    x,
    subsets,
    num.threads = 1,
    num.mads = 3,
    block = NULL,
    altexp.proportions = NULL,
    assay.type = "counts",
    output.prefix = NULL, 
    meta.name = "qc",
    flatten = TRUE
) {
    metrics <- computeRnaQcMetricsWithAltExps(x, subsets, altexp.proportions=altexp.proportions, num.threads=num.threads)
    thresholds <- scrapper::suggestRnaQcThresholds(metrics$main, block=block, num.mads=num.mads)
    keep <- scrapper::filterRnaQcMetrics(thresholds, metrics$main, block=block)

    df <- formatComputeRnaQcMetricsResult(metrics$main, flatten=flatten)
    df$keep <- keep
    colnames(df) <- paste0(output.prefix, colnames(df))
    colData(x) <- cbind(colData(x), df)

    if (!is.null(altexp.proportions)) {
        for (ae.name in names(metrics$altexp)) {
            ae.df <- formatComputeRnaQcMetricsResult(metrics$altexp[[ae.name]], flatten=flatten)
            colnames(ae.df) <- paste0(output.prefix, colnames(ae.df))
            colData(altExp(x, ae.name)) <- cbind(colData(altExp(x, ae.name)), ae.df)
        }
    }

    if (!is.null(meta.name)) {
        metadata(x)[[meta.name]] <- list(thresholds=thresholds)
    }

    x
}

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment altExp
#' @rdname quickRnaQc.se
computeRnaQcMetricsWithAltExps <- function(x, subsets, altexp.proportions, num.threads = 1, assay.type = "counts") {
    metrics <- scrapper::computeRnaQcMetrics(assay(x, assay.type, withDimnames=FALSE), subsets, num.threads=num.threads)

    # Adding more proportions from the alternative experiments, mostly for spike-ins.
    altexp.collected <- list()
    if (!is.null(altexp.proportions)) {
        altexp.proportions <- .sanitize_altexp_assays(altexp.proportions, all.altexps=altExpNames(x), default.assay.type=assay.type)
        for (ae.name in names(altexp.proportions)) {
            alt.assay <- assay(altExp(x, ae.name), altexp.proportions[[ae.name]])
            alt.metrics <- scrapper::computeRnaQcMetrics(alt.assay, subsets=list(), num.threads=num.threads)
            altexp.collected[[ae.name]] <- alt.metrics
            metrics$subsets[[ae.name]] <- alt.metrics$sum/(metrics$sum + alt.metrics$sum)
        }
    }

    list(main=metrics, altexp=altexp.collected)
}

#' @export
#' @rdname quickRnaQc.se
#' @importFrom S4Vectors DataFrame make_zero_col_DFrame
formatComputeRnaQcMetricsResult <- function(compute.res, flatten = TRUE) {
    df <- DataFrame(sum=compute.res$sum, detected=compute.res$detected)

    if (flatten) {
        for (sub in names(compute.res$subsets)) {
            df[[paste0(sub, ".proportion")]] <- compute.res$subsets[[sub]]
        }
    } else {
        tmp <- make_zero_col_DFrame(nrow=nrow(df))
        for (sub in names(compute.res$subsets)) {
            tmp[[sub]] <- compute.res$subsets[[sub]]
        }
        df[["proportion"]] <- tmp
    }

    df
}
