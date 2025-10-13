# We keep arguments that are explicitly named in 'more.args' for back-compatibility,
# in case some of these are promoted into top-level arguments at a later date.
.collapse_args <- function(named.args, more.args) {
    c(named.args[!(names(named.args) %in% names(more.args))], more.args)
}
