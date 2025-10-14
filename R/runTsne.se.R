#' t-stochastic neighbor embedding
#'
#' Generate a t-SNE visualization from an existing embedding.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param perplexity,num.threads Arguments to pass to \code{\link[scrapper]{runTsne}}.
#' @param more.tsne.args Named list of further arguments to pass to \code{\link[scrapper]{runTsne}}.
#' @param reddim.type Integer or string specifying the existing embedding in the \code{\link[SingleCellExperiment]{reducedDim}} of \code{x}.
#' @param output.name String containing the name of the output \code{\link[SingleCellExperiment]{reducedDim}}. 
#'
#' @return \code{x} is returned with the t-SNE coordinates stored in the \code{reducedDim}.
#'
#' @author Aaron Lun
#' @examples
#' sce <- getTestRnaData.se("pca")
#' # Using fewer iterations for a faster-running example.
#' sce <- runTsne.se(sce, more.tsne.args=list(max.iterations=50))
#' head(reducedDim(sce, "TSNE"))
#'  
#' @export
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
runTsne.se <- function(
    x,
    perplexity = 30,
    num.threads = 1,
    more.tsne.args = list(),
    reddim.type = "PCA", 
    output.name = "TSNE"
) {
    out <- do.call(
        scrapper::runTsne,
        c(
            list(t(reducedDim(x, reddim.type))),
            .collapse_args(
                list(perplexity=perplexity, num.threads=num.threads),
                more.tsne.args
            )
        )
    )

    reducedDim(x, output.name) <- out
    x
}
