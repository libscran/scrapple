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
#' @seealso
#' \code{\link[scrapper]{runTsne}} from the \pkg{scrapper} package.
#'  
#' @export
runTsne.se <- function(
    x,
    perplexity = 30,
    num.threads = 1,
    more.tsne.args = list(),
    reddim.type = "PCA", 
    output.name = "TSNE"
) {
    out <- .call(
        scrapper::runTsne,
        list(.get_transposed_reddim(x, reddim.type)),
        list(perplexity=perplexity, num.threads=num.threads),
        more.tsne.args
    )

    .add_tsne_results(x, output.name, out)
}

#' @importFrom SingleCellExperiment reducedDim<-
.add_tsne_results <- function(x, output.name, res) {
    reducedDim(x, output.name) <- res
    x
}
