#' @import magrittr


### Utility functions and macros.

get_single <- function (lst) {
    lst %T>% { stopifnot(length(.) == 1) } %>% .[[1]]
}

is_just_string <- function (x) {
    is.vector(x, mode = 'character') && (length(x) == 1)
}

## get names and perform some extra checks
get_names <- function (x, unique = TRUE) {
    names(x) %T>%
        { stopifnot(!is.null(.) || !any(. == '')) } %T>%
        { stopifnot(!unique || (anyDuplicated(.) == 0)) }
}

replace_colnames <- function (df, nm) {
    stopifnot(anyDuplicated(nm, incomparables = NA) == 0)
    nonna <- !is.na(nm)
    df[,nonna] %>% set_colnames(nm[nonna])
}

msg <- function (...) {
    message(sprintf(...))
}
