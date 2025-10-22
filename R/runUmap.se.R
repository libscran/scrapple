#' Uniform manifold approximation and projection
#'
#' Generate a UMAP visualization from an existing embedding.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param num.dim,num.neighbors,min.dist,num.threads Arguments to pass to \code{\link[scrapper]{runUmap}}.
#' @param more.umap.args Named list of further arguments to pass to \code{\link[scrapper]{runUmap}}.
#' @param reddim.type Integer or string specifying the existing embedding in the \code{\link[SingleCellExperiment]{reducedDim}} of \code{x}.
#' @param output.name String containing the name of the output \code{\link[SingleCellExperiment]{reducedDim}}. 
#'
#' @return \code{x} is returned with the UMAP coordinates stored in the \code{reducedDim}.
#'
#' @author Aaron Lun
#' @examples
#' sce <- getTestRnaData.se("pca")
#' # Using fewer epochs for a faster-running example.
#' sce <- runUmap.se(sce, more.umap.args=list(num.epochs=50))
#' head(reducedDim(sce, "UMAP"))
#'
#' @seealso
#' \code{\link[scrapper]{runUmap}} from the \pkg{scrapper} package.
#'  
#' @export
runUmap.se <- function(
    x,
    num.dim = 2,
    min.dist = 0.1,
    num.neighbors = 15,
    num.threads = 1,
    more.umap.args = list(),
    reddim.type = "PCA", 
    output.name = "UMAP"
) {
    res <- .call(
        scrapper::runUmap,
        list(.get_transposed_reddim(x, reddim.type)),
        list(num.dim=num.dim, min.dist=min.dist, num.neighbors=num.neighbors, num.threads=num.threads),
        more.umap.args
    )

    .add_umap_results(x, output.name, res)
}

#' @importFrom SingleCellExperiment reducedDim<-
.add_umap_results <- function(x, output.name, res) {
    reducedDim(x, output.name) <- res
    x
}
