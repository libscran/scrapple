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
