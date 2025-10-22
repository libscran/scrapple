#' Aggregate expression across gene sets
#'
#' Aggregate expression values across sets of genes for each cell.
#'
#' @param x A \link[SummarizedExperiment]{SummarizedExperiment} object or one of its subclasses.
#' Rows correspond to genes and columns correspond to cells.
#' @param sets,num.threads Arguments to pass to \code{\link[scrapper]{aggregateAcrossCells}}.
#' @param more.aggr.args Named list of additional arguments to pass to \code{\link[scrapper]{aggregateAcrossCells}}.
#' @param assay.type Integer or string specifying the assay of \code{x} to be aggregated.
#' @param output.name String specifying the assay name in the output object.
#' Defaults to \code{assay.type} if it is a string, otherwise \code{"aggregated"}.
#'
#' @details
#' \code{sets} may also be a \link[S4Vectors]{List} subclass,
#' in which case the \code{\link[S4Vectors]{mcols}} are used as the \code{\link[SummarizedExperiment]{rowData}} of the output object.
#' Weighted gene sets may be represented by a \link[IRanges]{DataFrameList} where each DataFrame contains two columns,
#' i.e., the gene identities and the associated weights.
#' 
#' @return A \link[SummarizedExperiment]{SummarizedExperiment} with number of rows equal to the number of gene sets.
#' The lone assay contains the aggregated values for each gene set for all cells. 
#' The \code{\link[SummarizedExperiment]{colData}} is the same as that of \code{x}.
#'
#' @author Aaron Lun
#' @examples
#' sce <- getTestRnaData.se("norm")
#'
#' library(org.Mm.eg.db)
#' some.sets <- select(
#'     org.Mm.eg.db,
#'     keytype="GO",
#'     keys=c(
#'         "GO:0048709", # oligodendrocyte differentiation 
#'         "GO:0048699", # neuron development
#'         "GO:0048143"  # astrocyte activation 
#'     ),
#'     columns="SYMBOL"
#' )
#' some.sets <- some.sets[some.sets$SYMBOL %in% rownames(sce),]
#' sets <- splitAsList(some.sets$SYMBOL, some.sets$GO)
#'
#' aggregated <- aggregateAcrossGenes.se(sce, sets)
#' aggregated
#' assay(aggregated)[,1:10]
#'
#' @export
#' @importFrom methods is
#' @importFrom S4Vectors mcols
#' @importClassesFrom S4Vectors List
#' @importClassesFrom IRanges DataFrameList
#' @importFrom SummarizedExperiment rowData<- colData SummarizedExperiment
aggregateAcrossGenes.se <- function(
    x,
    sets,
    num.threads = 1,
    more.aggr.args = list(),
    assay.type = "logcounts",
    output.name = NULL
) {
    rd <- NULL
    if (is(sets, "List")) {
        rd <- mcols(sets)
        if (!is.null(rd)) {
            rownames(rd) <- names(sets) # just in case
        }
        if (is(sets, "DataFrameList")) {
            sets <- as.list(sets)
            sets <- lapply(sets, as.list)
        } else {
            sets <- as.list(sets)
        }
    }

    vecs <- .call(
        scrapper::aggregateAcrossGenes,
        list(assay(x, assay.type)),
        list(sets=sets, num.threads=num.threads),
        more.aggr.args
    )
    if (length(vecs)) {
        mat <- do.call(rbind, vecs)
    } else {
        mat <- matrix(0, 0, ncol(x))
    }

    assays <- list(mat)
    if (!is.null(output.name)) {
        names(assays) <- output.name
    } else if (is.character(assay.type)) {
        names(assays) <- assay.type
    } else {
        names(assays) <- "aggregated"
    }

    out <- SummarizedExperiment(assays, colData=colData(x))
    if (!is.null(rd)) {
        rowData(out) <- rd
    }
    out
}
