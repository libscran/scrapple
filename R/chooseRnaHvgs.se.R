#' Choose highly variable genes from residuals
#'
#' Model the mean-variance relationship across genes and choose highly variable genes (HVGs) based on the residuals of the fitted trend.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param block,block.weight.policy,variable.block.weight,mean.filter,min.mean,transform,span,use.min.width,min.width,min.window.count,num.threads
#' Arguments passed to \code{\link[scrapper]{modelGeneVariances}}.
#' @param top,keep.ties,bound Arguments passed to \code{\link[scrapper]{chooseHighlyVariableGenes}}.
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
#' sce <- chooseRnaHvgs.se(sce, use.min.width=TRUE)
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
    block.weight.policy = "variable",
    variable.block.weight = c(0, 1000),
    mean.filter = TRUE,
    min.mean = 0.1,
    transform = TRUE,
    span = 0.3,
    use.min.width = FALSE,
    min.width = 1,
    min.window.count = 200,
    num.threads = 1,
    top = 4000,
    keep.ties = TRUE,
    bound = 0,
    assay.type = "logcounts",
    prefix = NULL,
    include.per.block = FALSE
) {
    info <- scrapper::modelGeneVariances(
        assay(x, assay.type),
        block=block,
        block.weight.policy=block.weight.policy,
        variable.block.weight=variable.block.weight,
        mean.filter=mean.filter,
        min.mean=min.mean,
        transform=transform,
        span=span,
        use.min.width=use.min.width,
        min.width=min.width,
        min.window.count=min.window.count,
        num.threads=num.threads
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

    is.hvg <- logical(nrow(x))
    hvg.index <- scrapper::chooseHighlyVariableGenes(info$statistics$residuals, top=top, larger=TRUE, keep.ties=keep.ties, bound=bound)
    is.hvg[hvg.index] <- TRUE
    rowData(x)[[paste0(prefix, "hvg")]] <- is.hvg

    x
}
