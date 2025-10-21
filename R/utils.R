# We keep arguments that are explicitly named in 'more.args' for back-compatibility,
# in case some of these are promoted into top-level arguments at a later date.
.collapse_args <- function(named.args, more.args) {
    c(named.args[!(names(named.args) %in% names(more.args))], more.args)
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
