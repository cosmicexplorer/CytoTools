#' @import magrittr


### Utility functions and macros.

get_single <- function (lst) {
    lst %T>% { stopifnot(length(.) == 1) } %>% .[[1]]
}

is_just_string <- function (x) {
    is.vector(x, mode = 'character') && (length(x) == 1)
}

is_counting_num <- function (x, allow_zero = FALSE) {
    is.vector(x, mode = 'integer') &&
        (length(x) == 1) &&
        ((x > 0) || (allow_zero && (x == 0)))
}

is_nonneg_scalar <- function (x, allow_zero = TRUE) {
    is.vector(x, mode = 'double') &&
        (length(x) == 1) &&
        ((x > 0) || (allow_zero && (x == 0)))
}

#' @title Rotate a Matrix Counter-Clockwise
#'
#' @description `rotate_matrix_ccw` performs a counter-clockwise rotation of the
#'     input matrix.
#'
#' @param mat Any matrix.
#'
#' @return Return an n x m matrix obtained by rotating `mat` counter-clockwise,
#'     if `mat` is an m x n matrix.
#'
#' @examples
#' mat <- matrix(1:8, ncol = 2)
#' mat
#' ##      [,1] [,2]
#' ## [1,]    1    5
#' ## [2,]    2    6
#' ## [3,]    3    7
#' ## [4,]    4    8
#' rotate_matrix_ccw(mat)
#' ##      [,1] [,2] [,3] [,4]
#' ## [1,]    5    6    7    8
#' ## [2,]    1    2    3    4
#'
#' @export
#'
rotate_matrix_ccw <- function (mat) { mat[,ncol(mat):1] %>% t }

#' @title Execute Expression and Measure Execution Time
#'
#' @description `timed_execute` executes an expression and measures the time
#'     taken, returning both.
#'
#' @param expr R expression to execute.
#' @param digits integer, number of digits to round time to.
#' @param time_type integer vector indexing into the result of [system.time()],
#'     defaulting to "wall clock time".
#' @param env environment to execute `expr` in, defaulting to the calling
#'     environment.
#'
#' @return Named list with `time` rounded to `digits` in seconds, `value` as the
#'     result of evaluating `expr`.
#'
#' @seealso [system.time()] is called to time the execution of `expr`.
#'
#' @examples
#' b <- 3
#' timed_execute({Sys.sleep(.5); b + 2})
#' ## $time
#' ## elapsed
#' ##     0.5
#' ##
#' ## $value
#' ## [1] 5
#'
#' @export
#'
timed_execute <- function (expr,
                           digits = 3, time_type = 3,
                           env = parent.frame()) {
    ex <- substitute(expr)
    time <- system.time(gcFirst = FALSE, expr = {
        value <- eval(ex, env)
    }) %>% { round(.[time_type], digits = digits) }
    list(time = time, value = value)
}

## check if names exist and are defined
check_names <- function (x, nonempty = TRUE, unique = TRUE) {
    is.vector(x, mode = 'character') &&
        (!nonempty || (length(x) > 0)) &&
        (!any(x == '')) &&
        (!unique || (anyDuplicated(x) == 0))
}

compare_names <- function (x, y, ...) {
    (check_names(x, ...) && check_names(y, ...)) &&
        (length(x) == length(y)) &&
        (x == y)
}

get_names <- function (x, ...) {
    names(x) %T>% { stopifnot(check_names(., ...)) }
}

## TODO: is this used anywhere?
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
            if (is.null(desc)) { "<no description>" } else { desc },
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
    stopifnot(is.list(a), is.list(b))
    if (length(a) == 0) { return(b) }
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

## TODO: is this used anywhere?
rotate_matrix_ccw <- function (m) { m[,ncol(m):1] %>% t }

safe_int_seq <- function (from, to) {
    stopifnot(is_counting_num(from),
              is_counting_num(to, allow_zero = TRUE))
    stopifnot(to >= (from - 1))
    if (to == from - 1) {
        integer(0)
    } else {
        from:to
    }
}
