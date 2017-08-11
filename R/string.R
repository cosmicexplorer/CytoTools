### String manipulation.

#' @import magrittr


#' @title Split Elements of a Character Vector By Multiple Elements
#'
#' @description `multi_str_split` splits the elements of `strs` into substrings
#'     between instances of any element of `splits`.
#'
#' @param strs character vector, each element of which is to be split.
#' @param splits character vector of fixed strings or PCRE regular expressions
#'     to use for splitting. If `splits` has length 0, `strs` is split
#'     into single characters. All elements of `splits` are used to split
#'     every element of `strs`!
#' @param fixed logical. If `TRUE` match `splits` exactly, otherwise
#'     use PCRE regular expressions.
#'
#' @return A list of the same length as `strs`, where each element is a
#'     character vector containing the splits of the corresponding element of
#'     `strs` according to `splits`.
#'
#' @seealso [strsplit()] is used to perform the splitting.
#'
#' @examples
#' multi_str_split(c("a_b-c d", "a-b_c d", "123"), c("-", "_", " ", "2"),
#'     fixed = T)
#' ## [[1]]
#' ## [1] "a" "b" "c" "d"
#' ##
#' ## [[2]]
#' ## [1] "a" "b" "c" "d"
#' ##
#' ## [[3]]
#' ## [1] "1" "3"
#'
#' @export
#'
multi_str_split <- function (strs, splits, fixed) {
    stopifnot(is.vector(strs, 'character'),
              is.vector(splits, 'character'),
              is.vector(fixed, 'logical'),
              length(fixed) == 1)
    if (length(splits) == 0) { return(strsplit(strs, character(0))) }
    do_split <- if (fixed) {
                    function (s, spl) strsplit(s, spl, fixed = TRUE)
                } else {
                    function (s, spl) strsplit(s, spl, perl = TRUE)
                }
    Reduce(x = splits, init = strs, f = function (s_list, spl) {
        lapply(s_list, (. %>% do_split(spl) %>% Reduce(f = c)))
    })
}

sort_component_helper <- function (splits, indices, orders) {
    n <- length(indices)
    if (n == 0) {
        return(indices)
    }
    empty_p <- lapply(indices, function (i) length(splits[[i]]) == 0) %>% unlist
    if (all(empty_p)) {
        return(indices)
    }
    nonempty_inds <- indices[!empty_p]
    next_orders <- orders[-1]
    next_splits <- lapply(splits, function (cur) cur[-1])
    cur_splits <- lapply(indices, function (i) splits[[i]][1]) %>% unlist
    recognized <- if (length(orders) > 0) { orders[[1]] } else { character() }
    stopifnot(!any(duplicated(recognized)))
    lvls <- cur_splits[!(cur_splits %in% recognized)] %>% sort %>% unique %>% {
        as.character(c(recognized, .))
    }
    nonempty_reduced <- Reduce(
        x = lvls,
        init = list(result = list(),
                    remaining = nonempty_inds),
        f = function (cur, level) {
            if (length(cur$remaining) == 0) {
                return(cur)
            }
            matching <- lapply(cur$remaining, function (i) {
                splits[[i]][1] == level
            }) %>% unlist
            matched_sorted <- sort_component_helper(
                next_splits, cur$remaining[matching], next_orders)
            list(result = c(cur$result, matched_sorted),
                 remaining = cur$remaining[!matching])
        })
    stopifnot(length(nonempty_reduced$remaining) == 0)
    c(indices[empty_p], nonempty_reduced$result) %>% unlist
}


#' @title Sort strings by splitting and ordering.
#'
#' @description `sort_by_component` splits strings into lists of
#'     components, then orders them lexicographically by the value of each
#'     succeeding component.
#'
#' @param strs character vector of strings to be sorted.
#' @param split_by character vector used to split `strs` into
#'     components. This is interpreted as either fixed strings or PCRE regular
#'     expressions according to `fixed`.
#' @param orders list of character vectors indicating a lexicographic ordering
#'     for each component.
#' @param fixed logical. How to interpret `split_by`. If `TRUE`,
#'     `split_by` is treated as fixed strings, otherwise as PCRE regular
#'     expressions.
#' @param value logical. Whether to return the sorted strings or the permutation
#'     of the input.
#'
#' @return If `value = TRUE` return the sorted character vector,
#'     otherwise an integer vector which represents the permutation of
#'     `strs` which produces the sorted result (as in [order()]).
#'
#' @details [multi_str_split()] is used to split each string of
#'     `strs` into a vector of substrings using the elements of
#'     `split_by`. `split_by` is interpreted as either fixed strings
#'     or PCRE regular expressions according to `fixed`.
#'
#'     Splitting a string produces a vector of substrings which occur between
#'     each occurrence of each splitting string or regular expression. For
#'     example, splitting the string `"A:B:C"` by `":"` produces the
#'     vector `c("A", "B", "C")`, which we call its *components*.
#'
#'     Note that `split_by = ""` or `split_by = character(0)` splits
#'     each string of `strs` into their individual characters as
#'     components.
#'
#'     Each succeeding character vector in `orders` is used to sort the
#'     component at the corresponding index in a lexicographic order. This means
#'     that the character vector at index 1 of `orders` is used to sort the
#'     first component of each input string, index 2 sorts the second, and so
#'     on, and that the sorting of earlier components takes precedence over
#'     later components.
#'
#'     The character vectors in `orders` are used to sort each component
#'     `c` of each string `s` at a given index `i` as
#'     follows. The components which match *exactly* some element of the
#'     vector `v` at index `i` of `orders` are sorted in the
#'     order given by `v`, while the components which do not show up in
#'     `v` are sorted using [order()]. The components which
#'     *do* show up in `v` are then sorted before the components which
#'     do not.
#'
#'     Note that this function will throw an exception if `v`
#'     contains any duplicates. In addition, if `length(orders) == n` for
#'     some integer `n`, all components at indices greater than `n`
#'     are sorted purely by their relative string value, but still only compared
#'     to other components at the same index.
#'
#'     Lexicographical (sometimes called hierarchical) sorting means that if two
#'     strings are split into the exact same components up to the index
#'     `k`, then differ at index `k + 1`, then their relative order is
#'     determined by the ordering of component `k + 1` of the two
#'     strings. If there is no component at which strings differ in sorting,
#'     they have no relative order, and they may be sorted either way. This
#'     means that this function does not perform a *stable sort*.
#'
#'     If one string `s_1` is split into `k` components and another
#'     string `s_2` splits into more than `k` components, and the two
#'     strings have equal components up to index `k`, then `s_1` is
#'     sorted before `s_2`.
#'
#' @seealso [multi_str_split()] is used to split each input string. [order()] is
#'     an example of another function which returns a permutation of the input
#'     as an integer vector.
#'
#' @examples
#' data_files <- c("MB004_6m_panel2.fcs", "MB004_3wk_panel2.fcs",
#'                 "MB004_12wk_panel2.fcs", "MB004_pre_panel2.fcs",
#'                 "MB005_12wk_panel2.fcs")
#'
#' sort_by_component(
#'     data_files,
#'     split_by = c("_"),
#'     orders = list(c("MB005"), c("pre", "3wk", "12wk", "6m")))
#'
#' ## [1] "MB005_12wk_panel2.fcs" "MB004_pre_panel2.fcs"
#' ## [3] "MB004_3wk_panel2.fcs"  "MB004_12wk_panel2.fcs"
#' ## [5] "MB004_6m_panel2.fcs"
#'
#' @export
#'
sort_by_component <- function (strs, split_by, orders = list(),
                               fixed = TRUE, value = TRUE) {
    splits <- multi_str_split(strs, split_by, fixed = fixed)
    indices <- sort_component_helper(splits, 1:length(splits), orders)
    if (value) { strs[indices] } else { indices }
}



sort_by_substrs <- function (strs) {
    
}
