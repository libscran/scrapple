# We keep arguments that are explicitly named in 'more.args' for back-compatibility,
# in case some of these are promoted into top-level arguments at a later date.
.call <- function(FUN, essential.args, named.args, more.args) {
    all.args <- c(
        essential.args,
        named.args[!(names(named.args) %in% names(more.args))],
        more.args
    )
    do.call(FUN, all.args)
}

.sanitize_altexp_assays <- function(altexps, all.altexps, default.assay.type) {
    if (!is.null(names(altexps))) {
        altexps[!duplicated(names(altexps))]
    } else {
        altexps <- unique(altexps)
        output <- rep(default.assay.type, length(altexps))
        if (is.numeric(altexps)) {
            names(output) <- all.altexps[altexps]
        } else {
            names(output) <- altexps
        }
        output
    }
}

#' @importFrom Matrix t
#' @importFrom DelayedArray DelayedArray
#' @importFrom SingleCellExperiment reducedDim<-
.add_transposed_reddim <- function(x, name, mat, delayed) {
    if (delayed) {
        mat <- DelayedArray(mat)
    }
    reducedDim(x, name) <- t(mat)
    x
}

#' @importFrom Matrix t
#' @importFrom methods is
#' @importClassesFrom DelayedArray DelayedArray
#' @importFrom SingleCellExperiment altExp reducedDim
.get_transposed_reddim <- function(x, name) {
    if (is.null(names(name))) {
        mat <- reducedDim(x, name, withDimnames=FALSE)
    } else {
        mat <- reducedDim(altExp(x, names(name), withDimnames=FALSE), name, withDimnames=FALSE)
    }

    mat <- t(mat)

    if (is.matrix(mat)) {
        colnames(mat) <- colnames(x)
        return(mat)
    } else if (is(mat, "DelayedArray")) {
        # Possibly a no-op if .add_transposed_reddim was set with delayed=TRUE.
        if (is.matrix(mat@seed)) {
            out <- mat@seed
            colnames(out) <- colnames(x)
            return(out)
        }
    }

    colnames(mat) <- colnames(x)
    return(as.matrix(mat))
}
