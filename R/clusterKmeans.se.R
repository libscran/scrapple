#' k-means clustering of cells
#'
#' Perform k-means clustering on an existing low-dimensional embedding.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param k Arguments to be passed to \code{\link[scrapper]{runKmeans}}.
#' @param more.kmeans.args Named list of further arguments to be passed to \code{\link[scrapper]{runKmeans}}.
#' @param reddim.type Integer or string specifying the existing embedding in the \code{\link[SingleCellExperiment]{reducedDim}} of \code{x}.
#' @param output.name String containing the name of the \code{\link[SummarizedExperiment]{colData}} column in which to store the cluster assignments.
#' @param meta.name String containing the name of the \code{\link[SummarizedExperiment]{metadata}} entry in which to store extra clustering output.
#' If \code{NULL}, no extra clustering output is stored. 
#'
#' @return \code{x} is returned with the cluster assignment for each cell stored in the \code{colData}.
#' Additional clustering output is stored in the \code{metadata}.
#'
#' @author Aaron Lun
#' @examples
#' sce <- getTestRnaData.se("pca")
#' sce <- clusterKmeans.se(sce, k=10)
#' table(sce$clusters)
#'
#' @seealso
#' \code{\link[scrapper]{clusterKmeans}} from the \pkg{scrapper} package.
#' 
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<-
clusterKmeans.se <- function(x, k, more.kmeans.args = list(), reddim.type = "PCA", output.name = "clusters", meta.name = NULL) {
    clout <- .call(
        scrapper::clusterKmeans,
        list(.get_transposed_reddim(x, reddim.type)),
        list(k=k),
        more.kmeans.args
    )

    colData(x)[[output.name]] <- clout$clusters
    if (!is.null(meta.name)) {
        clout$clusters <- NULL
        metadata(x)[[meta.name]] <- clout
    }

    x
}
