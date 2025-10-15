#' Quick quality control for RNA data
#'
#' Quickly compute quality control (QC) metrics, thresholds and filters from RNA data in a \link[SummarizedExperiment]{SummarizedExperiment}.
#' 
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param subsets,num.threads Arguments passed to \code{\link[scrapper]{computeRnaQcMetrics}}.
#' @param num.mads,block Arguments passed to \code{\link[scrapper]{suggestRnaQcThresholds}}.
#' @param altexp.proportions Character vector containing the names of alternative experiments for which to compute proportions relative to the RNA total.
#' These proportions will be used for filtering in the same manner as the proportions computed from \code{subsets}.
#' @param assay.type Integer or string specifying the assay of \code{x} (and its alternative experiments) containing the RNA count matrix.
#' @param output.prefix String containing a prefix to add to the names of the \code{link[SummarizedExperiment]{colData}} columns containing the output statistics.
#' @param thresholds.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry containing the filtering thresholds.
#' If \code{NULL}, thresholds are not stored in the metadata.
#' @param flatten Logical scalar indicating whether to flatten the subset proportions into separate columns of the \code{link[SummarizedExperiment]{colData}}.
#' If \code{FALSE}, the subset proportions are stored in a nested \link[S4Vectors]{DataFrame}.
#' @param compute.res List returned by \code{\link[scrapper]{computeRnaQcMetrics}}.
#' 
#' @details
#' \code{altexp.proportions} is typically used to refer to alternative experiments holding spike-in data.
#' For each alternative experiment, the proportion is defined \deqn{X/(X+Y)} where \deqn{X} is the alternative experiment's total and \deqn{Y} is the RNA total.
#' The count matrix in each alternative experiment should be stored in the assay specified by \code{assay.type}.
#' If the vector itself is named, those names are used in the column names of the output proportions, rather than the names of the alternative experiments.
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
#' colData(altExp(sce, "ERCC"))[,c("sum", "detected")]
#'
#' @export
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SingleCellExperiment altExp altExp<-
quickRnaQc.se <- function( 
    x,
    subsets,
    num.threads = 1,
    num.mads = 3,
    block = NULL,
    altexp.proportions = NULL,
    assay.type = "counts",
    output.prefix = NULL, 
    thresholds.name = "thresholds",
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
        altexp.names <- .choose_altexp_names(altexp.proportions)
        for (i in seq_along(altexp.proportions)) {
            ap <- altexp.proportions[i]
            ae.name <- altexp.names[i]
            ae.df <- formatComputeRnaQcMetricsResult(metrics$altexp[[ae.name]], flatten=flatten)
            colnames(ae.df) <- paste0(output.prefix, colnames(ae.df))
            colData(altExp(x, ap)) <- cbind(colData(altExp(x, ap)), ae.df)
        }
    }

    if (!is.null(thresholds.name)) {
        metadata(x)[[thresholds.name]] <- thresholds
    }

    x
}

#' @export
#' @importFrom SummarizedExperiment assay
#' @importFrom SingleCellExperiment altExp
computeRnaQcMetricsWithAltExps <- function(x, subsets, altexp.proportions, num.threads = 1, assay.type = "counts") {
    metrics <- scrapper::computeRnaQcMetrics(assay(x, assay.type, withDimnames=FALSE), subsets, num.threads=num.threads)

    # Adding more proportions from the alternative experiments, mostly for spike-ins.
    altexp.collected <- list()
    if (!is.null(altexp.proportions)) {
        altexp.names <- .choose_altexp_names(altexp.proportions)
        for (i in seq_along(altexp.proportions)) {
            alt.assay <- assay(altExp(x, altexp.proportions[i]), withDimnames=FALSE)
            alt.metrics <- scrapper::computeRnaQcMetrics(alt.assay, subsets=list(), num.threads=num.threads)
            ae.name <- altexp.names[i]
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

.choose_altexp_names <- function(altexp.proportions) {
    altexp.names <- names(altexp.proportions)
    if (is.null(altexp.names)) {
        altexp.proportions
    } else {
        altexp.names
    }
}
