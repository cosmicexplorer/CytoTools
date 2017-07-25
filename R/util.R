#' @import magrittr


### Utility functions and macros.

get_single <- function (lst) {
    lst %T>% { stopifnot(length(.) == 1) } %>% .[[1]]
}

## TODO: generate a switch by macro
## gen_switch <- function (...) {
##     do.call('switch', )
## }

is_just_string <- function (x) {
    is.vector(x, mode = 'character') && (length(x) == 1)
}

## get names and perform some extra checks
get_names <- function (x, unique = TRUE) {
    names(x) %T>%
        { stopifnot(!is.null(.) || !any(. == '')) } %T>%
        { stopifnot(!unique || !any(duplicated(.))) }
}
