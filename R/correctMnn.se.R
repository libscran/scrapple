#' MNN correction
#'
#' Correct batch effects from an existing embedding with mutual nearest neighbors (MNNs).
#'
#' @param x A \link[SingleCellExperiment]{SingleCellExperiment} object or one of its subclasses.
#' Rows correspond to genomic features and columns correspond to cells.
#' @param block Arguments passed to \code{\link[scrapper]{correctMnn}}.
#' @param more.mnn.args Named list of additional arguments to pass to \code{\link[scrapper]{correctMnn}}.
#' @param reddim.type String or integer specifying the \code{\link[SingleCellExperiment]{reducedDim}} entry on which to perform MNN correction.
#' @param output.name String containing the name of the \code{\link[SingleCellExperiment]{reducedDim}} entry in which to store the corrected embedding.
#'
#' @return \code{x} is returned with the corrected embedding stored as a \code{reducedDim} entry.
#' @author Aaron Lun
#'
#' @examples
#' sce <- getTestRnaData.se("pca")
#' # Treating the tissue of origin as the batch.
#' sce <- correctMnn.se(sce, sce$tissue)
#' reducedDimNames(sce)
#'
#' @seealso
#' \code{\link[scrapper]{correctMnn}} from the \pkg{scrapper} package.
#'
#' @export
#' @importFrom SingleCellExperiment reducedDim reducedDim<-
correctMnn.se <- function(
    x,
    block,
    more.mnn.args = list(),
    reddim.type = "PCA", 
    output.name = "MNN"
) {
    out <- do.call(
        scrapper::correctMnn,
        c(
            list(t(reducedDim(x, reddim.type))),
            .collapse_args(
                list(block=block),
                more.mnn.args
            )
        )
    )

    reducedDim(x, output.name) <- t(out$corrected)
    x
}
