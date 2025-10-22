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
#' @param delayed.transpose Logical scalar indicating whether to delay the transposition when storing coordinates in the \code{\link[SingleCellExperiment]{reducedDims}}.
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
correctMnn.se <- function(
    x,
    block,
    more.mnn.args = list(),
    reddim.type = "PCA", 
    output.name = "MNN",
    delayed.transpose = FALSE
) {
    out <- .call(
        scrapper::correctMnn,
        list(.get_transposed_reddim(x, reddim.type)),
        list(block=block),
        more.mnn.args
    )

    .add_transposed_reddim(x, output.name, out$corrected, delayed.transpose)
}
