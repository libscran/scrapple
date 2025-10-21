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
#' @param num.neighbors,num.threads Arguments to pass to \code{\link[scrapper]{scaleByNeighbors}}.
#' @param more.scale.args Named list of additional arguments to pass to \code{\link[scrapper]{scaleByNeighbors}}.
#' @param output.name String containing the name of the \code{\link[SingleCellExperiment]{reducedDim}} entry in which to store the combined embeddings.
#' @param meta.name String containing the name of the \code{\link[S4Vectors]{metadata}} entry in which to store additional metrics.
#' If \code{NULL}, additional metrics are not stored.
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
#' @importFrom SingleCellExperiment altExp altExpNames reducedDim reducedDimNames mainExpName reducedDim<-
#' @importFrom S4Vectors metadata metadata<-
scaleByNeighbors.se <- function(
    x,
    altexp.reddims,
    main.reddims = "PCA", 
    num.neighbors = 20,
    num.threads = 1,
    more.scale.args = list(),
    output.name = "combined",
    meta.name = "combined"
) {
    all.embeddings <- list()
    main.reddims <- .sanitize_reddims(main.reddims, reducedDimNames(x))
    for (r in main.reddims) {
        all.embeddings <- append(all.embeddings, list(t(reducedDim(x, r))))
    }

    altexp.reddims <- altexp.reddims[!duplicated(names(altexp.reddims))]
    for (ae in names(altexp.reddims)) {
        ae.se <- altExp(x, ae)
        cur.reddim <- .sanitize_reddims(altexp.reddims[[ae]], reducedDimNames(ae.se))
        for (r in cur.reddim) {
            all.embeddings <- append(all.embeddings, list(t(reducedDim(ae.se, r))))
        }
        altexp.reddims[[ae]] <- cur.reddim
    }

    out <- do.call(
        scrapper::scaleByNeighbors,
        c(
            list(all.embeddings),
            .collapse_args(
                list(num.neighbors=num.neighbors, num.threads=num.threads),
                more.scale.args
            )
        )
    )

    reducedDim(x, output.name) <- t(out$combined)
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
