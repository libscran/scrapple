#' Graph-based clustering of cells
#'
#' Construct a shared-nearest neighbor (SNN) graph from an existing low-dimensional embedding,
#' and apply community detection algorithms to obtain clusters of cells.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param num.neighbors,num.threads Arguments to be passed to \code{\link[scrapper]{buildSnnGraph}}.
#' @param more.snn.args Named list of further arguments to be passed to \code{\link[scrapper]{buildSnnGraph}}.
#' @param method,resolution Arguments to be passed to \code{\link[scrapper]{clusterGraph}}.
#' For \code{resolution}, this is either passed to \code{multilevel.resolution} or \code{leiden.resolution} depending on \code{method}.
#' @param more.graph.args Named list of further arguments to be passed to \code{\link[scrapper]{clusterGraph}}.
#' @param reddim.type Integer or string specifying the existing embedding in the \code{\link[SingleCellExperiment]{reducedDim}} of \code{x}.
#' @param output.name String containing the name of the column of the \code{\link[SummarizedExperiment]{colData}} in which to store the cluster assignments.
#' @param meta.name String containing the name of the \code{\link[SummarizedExperiment]{metadata}} entry in which to store extra clustering output.
#' If \code{NULL}, no extra clustering output is stored. 
#' @param graph.name String containing the name of the \code{\link[SummarizedExperiment]{metadata}} entry in which to store the SNN graph.
#' If \code{NULL}, the SNN graph is not stored. 
#'
#' @return \code{x} is returned with the cluster assignment for each cell stored in the \code{colData}.
#' Additional clustering output is stored in the \code{metadata}.
#'
#' @author Aaron Lun
#' @examples
#' sce <- getTestRnaData.se("pca")
#' sce <- clusterGraph.se(sce)
#' table(sce$clusters)
#'
#' @seealso
#' \code{\link[scrapper]{buildSnnGraph}} and \code{\link[scrapper]{clusterGraph}}, from the \pkg{scrapper} package.
#' 
#' @export
#' @importFrom SingleCellExperiment reducedDim
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom SummarizedExperiment colData colData<-
clusterGraph.se <- function(
    x,
    num.neighbors = 10,
    num.threads = 1,
    more.snn.args = list(),
    method = "multilevel",
    resolution = NULL,
    more.graph.args = list(),
    reddim.type = "PCA",
    output.name = "clusters",
    meta.name = NULL,
    graph.name = NULL
) {
    graph.out <- do.call(
        scrapper::buildSnnGraph,
        c(
            list(t(reducedDim(x, reddim.type))),
            .collapse_args(
                list(num.neighbors=num.neighbors, num.threads=num.threads, as.pointer=is.null(graph.name)),
                more.snn.args
            )
        )
    )
    if (!is.null(graph.name)) {
        metadata(x)[[graph.name]] <- graph.out
    }

    res.args <- list()
    if (!is.null(resolution)) {
        res.args$multilevel.resolution <- resolution
        res.args$leiden.resolution <- resolution
    }

    clust.out <- do.call(
        scrapper::clusterGraph,
        c(
            list(graph.out),
            .collapse_args(res.args, more.graph.args)
        )
    )

    colData(x)[[output.name]] <- clust.out$membership
    if (!is.null(meta.name)) {
        clust.out$membership <- NULL
        metadata(x)[[meta.name]] <- clust.out
    }

    x
}
