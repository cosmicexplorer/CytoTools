### Routines for manipulating CyToF data and interacting with Cytobank.
## Written by Danny McClanahan, Irish Lab June 2017.
## <danieldmcclanahan@gmail.com>

#' @import magrittr
#' @import dplyr


### Read/write different representations of flow data.

#' @title List Files in a Directory
#'
#' @description `list_dir_files` lists all files in a directory without
#'     including other directories and using the full filenames.
#'
#' @param path string, the directory to read.
#'
#' @return Character vector containing paths to all files located at the
#'     specified `path`.
#'
#' @examples
#' ## list all files in the current directory
#' list_dir_files(".")
#'
#' @export
#'
list_dir_files <- function (path) {
    list.files(path = path, pattern = NULL,
               all.files = TRUE, full.names = TRUE, recursive = FALSE,
               include.dirs = FALSE, no.. = TRUE)
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
#' @description `read_cyto_file` reads a file containing CyToF data into a
#'     data frame.
#'
#' @param fname character vector of length one. The filename to read in.
#' @param allow_skip logical indicating whether to allow reading text files with
#'     one blank line at the top.
#'
#' @details Filenames can be binary FCS files or text files with headers. Files
#'     with extension ".fcs" will be read as FCS files with
#'     [flowCore::read.FCS()]. ".txt" files will be read as TSV, while
#'     ".csv" files will be read as CSV. Extensions are interpreted
#'     case-insensitively, but files with unrecognized extensions will trigger
#'     an exception.
#'
#'     If the file is a text file (.txt or .csv) and starts with a blank line,
#'     this function will recognize that and skip the initial blank line, unless
#'     `allow_skip = FALSE`. This handles a known quirk in many real-world
#'     datasets.
#'
#'     The [flowCore::exprs()] function is used when turning an FCS
#'     file into a data frame, which drops metadata.
#'
#' @return [read_cyto_file()] returns a data frame containing the
#'     content of the specified data file.
#'
#' @seealso [read.table()] is used to read TSV and CSV files.
#'     [flowCore::read.FCS()] and [flowCore::exprs()] are
#'     used to read FCS files.
#'
#' @examples
#' ## read the file "iPSCs.fcs" in the current directory into `fcs_df`
#' fcs_df <- read_cyto_file("./iPSCs.fcs")
#'
#' fcs_df[1:4,c(1:2,11:12)]
#' ##   Time Cell_length        140   ICOS-141
#' ## 1   45          23 -1.0545764 -0.3860820
#' ## 2  442          34 -0.2881902 -1.1049303
#' ## 3  480          23 -1.0798522 -0.4776470
#' ## 4  599          38 -0.4168313 -0.6641003
#'
#' @export
#'
read_cyto_file <- function (fname, allow_skip = TRUE) {
    ## TODO: consider having a cache for this function if files are reused a lot
    ext <- tools::file_ext(fname) %>% tolower
    switch(
        ext,
        fcs = read_fcs_cyto_frame(fname),
        csv = read_text_cyto_frame(fname, allow_skip, sep = ","),
        txt = read_text_cyto_frame(fname, allow_skip, sep = "\t"),
        stop(sprintf("unrecognized extension '%s' for file '%s'",
                     ext, fname)))
}

#' @title Scan a Directory for CytoF Files, Sort Them, and Read Them Into Memory
#'
#' @description `read_files_sorted` reads all files with the given extensions in
#'     the current directory into a named list of data frames.
#'
#' @inheritParams list_dir_files
#' @param extensions character vector containing file extensions to process.
#' @inheritParams sort_by_component
#' @inheritParams read_cyto_file
#'
#' @details The directory specified by `path` is scanned for files with any of
#'     the given `extensions`. The elements of `split_by` are interpreted as
#'     *fixed strings* to split the resulting filenames, and `orders` are used
#'     to lexicographically sort them with [sort_by_component()]. The resulting
#'     sorted character vector of filenames is then used as the [names()] for a
#'     named list, which is then populated with each file's contents using
#'     [read_cyto_file()].
#'
#' @return Named list of data frames, where each name is the filename from which
#'     the value at that index is read.
#'
#' @seealso [list_dir_files()] is used to scan the directory,
#'     [sort_by_component()] sorts file paths, and [read_cyto_file()] is called
#'     to interpret each file's contents into a data frame.
#'
#' @examples
#' sorted_cytof_data <- read_files_sorted(
#'    path = ".",
#'    extensions = c("fcs"),
#'    split_by = c("_", " "),
#'    orders = list(c(), c("pre", "3wk", "12wk", "6m")))
#'
#' length(sorted_cytof_data)
#' ## [1] 141
#' is.data.frame(sorted_cytof_data[[1]])
#' ## [1] TRUE
#'
#' lapply(sorted_cytof_data[1:2], function (fcs_df) fcs_df[1:4,c(1:2,5:6)])
#' ## $`./29_RCCPBMC_panel4_PD-1+ CD4.fcs`
#' ##    Time Cell_length        140  ICOS-141
#' ## 1 19648          24  0.2647168  1.667185
#' ## 2 20933          17  1.0358380  9.964522
#' ## 3 27051          24 -0.8677309 17.684467
#' ## 4 28387          37  1.2802999 14.215277
#' ##
#' ## $`./29_RCCPBMC_panel4_PD-1+ CD8.fcs`
#' ##    Time Cell_length        140    ICOS-141
#' ## 1 29018          22 -0.9847549 -0.51542509
#' ## 2 30535          31 -0.5153459 -0.65859687
#' ## 3 33917          19 -1.0417140  3.72425652
#' ## 4 42929          27  1.7391516  0.03845136
#'
#' @export
#'
read_files_sorted <- function (path = ".",
                               extensions = c("fcs"),
                               split_by = c("_", " "),
                               orders = list(),
                               allow_skip = TRUE) {
    desired_files <- list_dir_files(path) %>%
        .[lapply(., tools::file_ext) %in% extensions]
    specified_order <- desired_files %>% lapply(basename) %>% unlist %>%
        sort_by_component(split_by = split_by, orders = orders,
                          fixed = TRUE, value = FALSE)
    sorted_files <- desired_files[specified_order]
    sorted_files %>%
        lapply((. %>% read_cyto_file(allow_skip = allow_skip))) %>%
        set_names(sorted_files)
}
