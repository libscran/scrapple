#' Choose highly variable genes from residuals
#'
#' Model the mean-variance relationship across genes and choose highly variable genes (HVGs) based on the residuals of the fitted trend.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param block,num.threads Arguments passed to \code{\link[scrapper]{modelGeneVariances}}.
#' @param more.var.args Named list of arguments to pass to \code{\link[scrapper]{modelGeneVariances}}.
#' @param top Arguments passed to \code{\link[scrapper]{chooseHighlyVariableGenes}}.
#' @param more.choose.args Named list of arguments to pass to \code{\link[scrapper]{chooseHighlyVariableGenes}}.
#' @param assay.type Integer or string specifying the assay of \code{x} containing the log-normalized expression matrix for the RNA data.
#' @param prefix String containing a prefix to add to the column names in the output \code{link[SummarizedExperiment]{rowData}}.
#' @param include.per.block Logical scalar indicating whether the per-block statistics should be stored in the output object.
#' Only relevant if \code{block} is specified.
#'
#' @return \code{x} is returned with the per-gene variance modelling statistics added to the \code{rowData}.
#' The \code{hvg} column indicates whether a gene was chosen as a HVG.
#' If \code{include.per.block=TRUE} and \code{block} is specified, the per-block statistics are stored as a nested DataFrame in the \code{per.block} column.
#'
#' @examples
#' sce <- getTestRnaData.se("norm")
#' sce <- chooseRnaHvgs.se(sce, more.var.args=list(use.min.width=TRUE))
#' summary(rowData(sce)$hvg)
#'
#' plot(rowData(sce)$means, rowData(sce)$variances, col=factor(rowData(sce)$hvg))
#' curve(approxfun(rowData(sce)$means, rowData(sce)$fitted)(x), col="dodgerblue", add=TRUE)
#' 
#' @export
#' @importFrom SummarizedExperiment assay rowData rowData<-
chooseRnaHvgs.se <- function(
    x, 
    block = NULL,
    num.threads = 1,
    more.var.args = list(),
    top = 4000,
    more.choose.args = list(),
    assay.type = "logcounts",
    prefix = NULL,
    include.per.block = FALSE
) {
    info <- do.call(
        scrapper::modelGeneVariances,
        c(
            list(assay(x, assay.type)),
            .collapse_args(
                list(block=block, num.threads=num.threads),
                more.var.args
            )
        )
    )

    colnames(info$statistics) <- paste0(prefix, colnames(info$statistics))
    rowData(x) <- cbind(rowData(x), info$statistics)

    if (include.per.block && !is.null(info$per.block)) {
        tmp <- S4Vectors::make_zero_col_DFrame(nrow=nrow(x))
        for (n in names(info$per.block)) {
            tmp[[sub]] <- DataFrame(info$per.block[[n]])
        }
        rowData(x)[[paste0(prefix, "per.block")]] <- tmp
    }

    hvg.index <- do.call(
        scrapper::chooseHighlyVariableGenes,
        c(
            list(info$statistics$residuals),
            .collapse_args(
                list(top=top, larger=TRUE),
                more.choose.args
            )
        )
    )

    is.hvg <- logical(nrow(x))
    is.hvg[hvg.index] <- TRUE
    rowData(x)[[paste0(prefix, "hvg")]] <- is.hvg

    x
}
