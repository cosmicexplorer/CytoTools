### Routines for manipulating CyToF data and interacting with Cytobank.
## Written by Danny McClanahan, Irish Lab June 2017.
## <danieldmcclanahan@gmail.com>

#' @import magrittr
#' @import dplyr


### Read/write different representations of flow data.

#' @title ?
#'
#' @description ?
#'
#' @param path ?
#'
#' @return ?
#'
#' @export
#'
fcs_file_paths <- function (path = ".", pattern = "\\.fcs$") {
    list.files(path = path, pattern = pattern,
               ignore.case = TRUE, all.files = TRUE, full.names = TRUE,
               recursive = FALSE, no.. = TRUE)
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
#' @description \code{sort_by_component} splits strings into lists of
#'     components, then orders them lexicographically by the value of each
#'     succeeding component.
#'
#' @param strs character vector of strings to be sorted.
#' @param split_by character vector used to split \code{strs} into
#'     components. This is interpreted as either fixed strings or PCRE regular
#'     expressions according to \code{fixed}.
#' @param orders list of character vectors indicating a lexicographic ordering
#'     for each component.
#' @param fixed logical. How to interpret \code{split_by}. If \code{TRUE},
#'     \code{split_by} is treated as fixed strings, otherwise as PCRE regular
#'     expressions.
#' @param value logical. Whether to return the sorted strings or the permutation
#'     of the input.
#'
#' @return If \code{value = TRUE} return the sorted character vector,
#'     otherwise an integer vector which represents the permutation of
#'     \code{strs} which produces the sorted result (as in \code{\link{order}}).
#'
#' @details \code{\link{strsplit}} is used to split each string of \code{strs}
#'     into a vector of substrings using the elements of
#'     \code{split_by}. \code{split_by} is interpreted as either fixed strings
#'     or PCRE regular expressions according to \code{fixed}.
#'
#'     Splitting a string produces a vector of substrings which occur between
#'     each occurrence of each splitting string or regular expression. For
#'     example, splitting the string \code{"A:B:C"} by \code{":"} produces the
#'     vector \code{c("A", "B", "C")}, which we call its \emph{components}.
#'
#'     Note that \code{split_by = ""} or \code{split_by = character(0)} splits
#'     each string of \code{strs} into their individual characters as
#'     components.
#'
#'     Each succeeding character vector in \code{orders} is used to sort the
#'     component at the corresponding index in a lexicographic order. This means
#'     that the character vector at index 1 of \code{orders} is used to sort the
#'     first component of each input string, index 2 sorts the second, and so
#'     on, and that the sorting of earlier components takes precedence over
#'     later components.
#'
#'     The character vectors in \code{orders} are used to sort each component
#'     \code{c} of each string \code{s} at a given index \code{i} as
#'     follows. The components which match \emph{exactly} some element of the
#'     vector \code{v} at index \code{i} of \code{orders} are sorted in the
#'     order given by \code{v}, while the components which do not show up in
#'     \code{v} are sorted using \code{\link{order}}. The components which
#'     \emph{do} show up in \code{v} are then sorted before the components which
#'     do not.
#'
#'     Note that this function will throw an exception if \code{v}
#'     contains any duplicates. In addition, if \code{length(orders) == n} for
#'     some integer \code{n}, all components at indices greater than \code{n}
#'     are sorted purely by their relative string value, but still only compared
#'     to other components at the same index.
#'
#'     Lexicographical (sometimes called hierarchical) sorting means that if two
#'     strings are split into the exact same components up to the index
#'     \code{k}, then differ at index \code{k + 1}, then their relative order is
#'     determined by the ordering of component \code{k + 1} of the two
#'     strings. If there is no component at which strings differ in sorting,
#'     they have no relative order, and they may be sorted either way. This
#'     means that this function does not perform a \emph{stable sort}.
#'
#'     If one string \code{s_1} is split into \code{k} components and another
#'     string \code{s_2} splits into more than \code{k} components, and the two
#'     strings have equal components up to index \code{k}, then \code{s_1} is
#'     sorted before \code{s_2}.
#'
#' @seealso \code{\link{strsplit}} is used to split each input
#'     string. \code{\link{order}} is a function which returns a permutation of
#'     the input as an integer vector.
#'
#' @examples
#' data_files <- c("MB004_6m_panel2.fcs", "MB004_3wk_panel2.fcs",
#'                 "MB004_12wk_panel2.fcs", "MB004_pre_panel2.fcs",
#'                 "MB005_12wk_panel2.fcs")
#' sort_by_component(data_files,
#'                   split_by = c("_"),
#'                   orders = list(c("MB005"), c("pre", "3wk", "12wk", "6m")))
#' ## [1] "MB005_12wk_panel2.fcs" "MB004_pre_panel2.fcs"
#' ## [3] "MB004_3wk_panel2.fcs"  "MB004_12wk_panel2.fcs" "MB004_6m_panel2.fcs"
#'
#' @export
#'
sort_by_component <- function (strs, split_by, orders = list(),
                               fixed = TRUE, value = TRUE) {
    splits <- if (fixed) {
                  strsplit(strs, split_by, fixed = TRUE)
              } else {
                  strsplit(strs, split_by, perl = TRUE)
              }
    indices <- sort_component_helper(splits, 1:length(splits), orders)
    if (value) {
        strs[indices]
    } else {
        indices
    }
}

read_fcs_flowFrame <- function (fname) {
    flowCore::read.FCS(
        fname, transformation = NULL, truncate_max_range = FALSE) %T>%
        ## die if there's more than one dataset in the fcs file
        { stopifnot(.@description[["$NEXTDATA"]] == "0") }
}

read_fcs_cyto_frame <- function (fname) {
    read_fcs_flowFrame(fname) %>% flowCore::exprs(.) %>% as.data.frame
}

read_text_cyto_frame <- function (fname, allow_skip, ...) {
    tryCatch(
        read.table(fname, header = TRUE, ...),
        error = function (e) {
            if (!allow_skip) { stop(e) }
            read.table(fname, header = TRUE, skip = 1, ...)
        }
    )
}


#' @title Read a CyToF file.
#'
#' @description \code{read_cyto_file} reads a file containing CyToF data into a
#'     data frame.
#'
#' @param fname character vector of length one. The filename to read in.
#' @param rx_replace named char vector (which may be empty, or NULL), where
#'     names are regular expressions ("regexes") to match against CyToF marker
#'     names, and values are replacements.
#' @param allow_skip logical indicating whether to allow reading text files with
#'     one blank line at the top.
#'
#' @details Filenames can be binary FCS files or text files with headers. Files
#'     with extension ".fcs" will be read as FCS files with
#'     \code{\link{flowCore::read.FCS}}. ".txt" files will be read as TSV, while
#'     ".csv" files will be read as CSV. Extensions are interpreted
#'     case-insensitively, but files with unrecognized extensions will trigger
#'     an exception.
#'
#'     If the file is a text file (.txt or .csv) and starts with a blank line,
#'     this function will recognize that and skip the initial blank line, unless
#'     \code{allow_skip = FALSE}. This handles a known quirk in many real-world
#'     datasets.
#'
#' @return \code{\link{read_cyto_file}} returns a data frame containing the
#'     content of the specified data file. Only channel values are retained --
#'     any metadata or parameters are dropped.
#'
#' @seealso \code{\link{gsub}} for basic examples of regex replacement, while
#'     \code{\link{stringr::str_replace_all}} is the function called to perform
#'     these replacements.
#'
#'     \code{\link{flowCore::read.FCS}} is used to read FCS files, while
#'     \code{\link{read.table}} is used to read TSV and CSV files.
#'
#' @export
#'
read_cyto_file <- function (fname, rx_replace = NULL, allow_skip = TRUE) {
    ## TODO: does any kind of data cleaning make sense here? see ../README.md
    ## TODO: consider having a cache for this function if files are reused a lot
    ext <- tools::file_ext(fname) %>% tolower
    df <- switch(
        ext,
        fcs = read_fcs_cyto_frame(fname),
        csv = read_text_cyto_frame(fname, allow_skip, sep = ","),
        txt = read_text_cyto_frame(fname, allow_skip, sep = "\t"),
        stop(sprintf("unrecognized extension '%s' for file '%s'",
                     ext, fname)))
    cols <- colnames(df) %T>% { stopifnot(!any(duplicated(.))) }
    newcols <-
        if (is.null(rx_replace)) {
            cols
        } else {
            stopifnot(is.vector(rx_replace, 'character') &&
                      is.vector(get_names(rx_replace), 'character'))
            stringr::str_replace_all(cols, rx_replace)
        } %T>% {
            stopifnot((length(.) == length(cols)) &&
                      !any(duplicated(.)))
        }
    set_colnames(df, newcols)
}
