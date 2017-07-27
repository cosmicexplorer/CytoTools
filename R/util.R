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

#' @title ?
#'
#' @description ?
#'
#' @param strs ?
#' @param replace_spec ?
#'
#' @return ?
#'
#' @details ?
#'
#' @export
#'
replace_matches <- function (strs, replace_spec) {
    rx <- get_names(replace_spec)
    Reduce(x = 1:length(replace_spec), init = rep(strs), f = function (cur, i) {
        pat <- rx[[i]]
        rpl <- replace_spec[[i]]
        matches <- grepl(pat, cur, perl = T)
        cur[matches] <- if (is.character(rpl)) {
                            gsub(pat, rpl, cur[matches], perl = T)
                        } else if (is.function(rpl)) {
                            rpl(cur[matches])
                        } else {
                            rpl
                        }
        cur
    })
}

replace_colnames <- function (df, desc, nm) {
    if (anyDuplicated(nm, incomparables = NA) != 0) {
        stop(sprintf(paste(
            sep = "\n",
            "The channel names for population '%s' contain",
            "duplicates after normalizing. Either two channels",
            "are being merged into one by mistake, or an",
            "intercalator/non-phenotype channel is not being removed.",
            "Modify the `channel_name_ops` argument to fix this.",
            "The duplicates are: [%s]",
            "The whole list of normalized channel names is: [%s]"),
            desc,
            paste0(collapse = ", ", nm[duplicated(nm, incomparables = NA)]),
            paste0(collapse = ", ", nm)))
    }
    nonna <- !is.na(nm)
    df[,nonna] %>% set_colnames(nm[nonna])
}

msg <- function (...) {
    message(sprintf(...))
}

merge_named_lists <- function (a, b) {
    nma <- get_names(a)
    nmb <- get_names(b)
    ret <- rep(a)
    for (b_name in nmb) {
        if (!(b_name %in% nma)) {
            ret[[b_name]] <- b[[b_name]]
        }
    }
    ret
}
