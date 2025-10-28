#' Scale and combine multiple embeddings
#'
#' Scale embeddings for different modalities to equalize their intra-population variance,
#' and then combine them into a single embedding for downstream analysis.
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param altexp.reddims Named list of character or integer vectors.
#' Each entry is named after an alternative experiment.
#' Each vector contains the names/indices of the \code{\link[SingleCellExperiment]{reducedDim}} embeddings from that experiment to be combined.
#' @param main.reddims Character or integer vector specifying the names/indices of the \code{\link[SingleCellExperiment]{reducedDim}} entries from \code{x} to be combined.
#' @param num.neighbors,block,BNPARAM,num.threads Arguments to pass to \code{\link[scrapper]{scaleByNeighbors}}.
#' @param more.scale.args Named list of additional arguments to pass to \code{\link[scrapper]{scaleByNeighbors}}.
#' @param output.name String containing the name of the \code{\link[SingleCellExperiment]{reducedDim}} entry in which to store the combined embeddings.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry in which to store additional metrics.
#' If \code{NULL}, additional metrics are not stored.
#' @param delayed.transpose Logical scalar indicating whether to delay the transposition when storing coordinates in the \code{\link[SingleCellExperiment]{reducedDims}}.
#'
#' @return \code{x} is returned with the combined embeddings stored in its \code{rowData}.
#' The scaling factors for all embeddings are stored in the \code{metadata}.
#'
#' @author Aaron Lun
#'
#' @examples
#' sce <- getTestAdtData.se("pca")
#' sce <- scaleByNeighbors.se(sce, altexp.reddims=list(ADT="PCA"))
#' reducedDimNames(sce) 
#' metadata(sce)$combined
#'
#' @seealso
#' \code{\link{scaleByNeighbors}} from the \pkg{scrapper} package.
#'
#' @export
#' @importFrom SingleCellExperiment altExp altExpNames reducedDim reducedDimNames reducedDim<-
#' @importFrom S4Vectors metadata metadata<-
#' @importFrom BiocNeighbors AnnoyParam
scaleByNeighbors.se <- function(
    x,
    altexp.reddims,
    main.reddims = "PCA", 
    num.neighbors = 20,
    block = NULL,
    BNPARAM = AnnoyParam(),
    num.threads = 1,
    more.scale.args = list(),
    output.name = "combined",
    meta.name = "combined",
    delayed.transpose = FALSE
) {
    all.embeddings <- list()
    main.reddims <- .sanitize_reddims(main.reddims, reducedDimNames(x))
    for (r in main.reddims) {
        all.embeddings <- append(all.embeddings, list(.get_transposed_reddim(x, r)))
    }

    altexp.reddims <- altexp.reddims[!duplicated(names(altexp.reddims))]
    for (ae in names(altexp.reddims)) {
        ae.se <- altExp(x, ae)
        cur.reddim <- .sanitize_reddims(altexp.reddims[[ae]], reducedDimNames(ae.se))
        for (r in cur.reddim) {
            all.embeddings <- append(all.embeddings, list(.get_transposed_reddim(ae.se, r)))
        }
        altexp.reddims[[ae]] <- cur.reddim
    }

    out <- .call(
        scrapper::scaleByNeighbors,
        list(all.embeddings),
        list(num.neighbors=num.neighbors, block=block, num.threads=num.threads, BNPARAM=BNPARAM),
        more.scale.args
    )

    x <- .add_transposed_reddim(x, output.name, out$combined, delayed.transpose)
    if (!is.null(meta.name)) {
        out$combined <- NULL

        # Formatting it in the same manner as the arguments.
        counter <- 1L
        out$main.scaling <- numeric(0)
        for (r in main.reddims) {
            out$main.scaling[[r]] <- out$scaling[[counter]]
            counter <- counter + 1L
        }

        out$altexp.scaling <- list()
        for (ae in names(altexp.reddims)) {
            current <- numeric(0)
            for (r in altexp.reddims[[ae]]) {
                current[[r]] <- out$scaling[[counter]]
                counter <- counter + 1L
            }
            out$altexp.scaling[[ae]] <- current
        }

        out$scaling <- NULL
        metadata(x)[[meta.name]] <- out
    }

    x
}

.sanitize_reddims <- function(reddims, all.reddims) {
    reddims <- unique(reddims)
    if (is.integer(reddims)) {
        all.reddims[reddims]
    } else {
        reddims
    }
}
