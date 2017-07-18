### Routines for manipulating CyToF data and interacting with Cytobank.
## Written by Danny McClanahan, Irish Lab June 2017.
## <danieldmcclanahan@gmail.com>

#' @import magrittr
#' @import dplyr


### Libraries


### Utility functions and macros.

get_single <- function (lst) {
    lst %T>% { stopifnot(length(.) == 1) } %>% .[[1]]
}



### Read/write different representations of flow data.

read_fcs_flowFrame <- function (fname) {
    flowCore::read.FCS(
        fname, transformation = NULL, truncate_max_range = F) %T>%
        ## die if there's more than one dataset in the fcs file
        { stopifnot(.@description[["$NEXTDATA"]] == "0") }
}

read_fcs_cyto_frame <- function (fname) {
    read_fcs_flowFrame(fname) %>% flowCore::exprs(.) %>% as.data.frame
}

read_text_cyto_frame <- function (fname, allow_skip, ...) {
    tryCatch(
        read.table(fname, header = T, ...),
        error = function (e) {
            if (!allow_skip) { stop(e) }
            read.table(fname, header = T, skip = 1, ...)
        }
    )
}

## TODO: generate this instead of a list
## gen_switch <- function (...) {
##     do.call('switch', )
## }

is_just_string <- function (x) {
    is.vector(x, mode = 'character') && (length(x) == 1)
}

get_names <- function (x, unique = T) {
    names(x) %T>%
        { stopifnot(!is.null(.) || !any(. == '')) } %T>%
        { stopifnot(!unique || !any(duplicated(.))) }
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
#' @seealso \code{\link{gsub}} for basic examples of regex replacement, while
#'     \code{\link{stringr::str_replace_all}} is the function called to perform
#'     these replacements.
#'
#'     \code{\link{flowCore::read.FCS}} is used to read FCS files, while
#'     \code{\link{read.table}} is used to read TSV and CSV files.
#'
#' @rdname process_cyto_dataset
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

#' @title Read many CyToF files.
#'
#' @description \code{process_cyto_dataset} reads CyToF data files with
#'     \code{\link{read_cyto_file}}.
#'
#' @param fnames character vector of filenames to read in.
#' @param ... arguments to pass to \code{\link{read_cyto_file}}.
#'
#' @return A named list of data frames containing the content of each file in
#'     \code{fnames}. Names correspond to the filename which was read to produce
#'     each data frame.
#'
#' @export
#'
process_cyto_dataset <- function (fnames, ...) {
    lapply(fnames, function (file) read_cyto_file(file, ...)) %>%
        set_names(fnames)
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
                               fixed = T, value = T) {
    splits <- if (fixed) {
                  strsplit(strs, split_by, fixed = T)
              } else {
                  strsplit(strs, split_by, perl = T)
              }
    indices <- sort_component_helper(splits, 1:length(splits), orders)
    if (value) {
        strs[indices]
    } else {
        indices
    }
}



### Clean fcs data.

numeric_pat <- "^[[:digit:]]+$"
sne_pat <- "[sS][nN][eE]"

is_double_vec <- function (vec) {
    is.vector(vec, mode = "double")
}

is_all_integer <- function (vec) {
    ## this benchmarked as the fastest method
    all(vec == as.integer(vec))
}

## TODO: get a citation for asinh_transform and justify it. get alternatives and
## comparisons
#' @title Hyperbolic Arcsin Transform for Cytometry Data
#'
#' @description \code{asinh_transform} applies the standard transform applied to
#'     channels of CyToF or flow cytometry data when comparing channels
#'     directly.
#'
#' @param x vector to transform
#'
#' @return The transformed vector.
#'
#' @details Performs the \code{asinh(x / 5)} transformation to the input data.
#'
#' @references Garner, Duane L., et al. "Fluorometric assessments of
#'     mitochondrial function and viability in cryopreserved bovine
#'     spermatozoa." Biology of Reproduction 57.6 (1997): 1401-1406.
#'
#' @export
#'
asinh_transform <- function (x) { asinh(x / 5) }

## Return data frame which contains only marker data.
cyto_data_cols <- function (frame,
                            excl_pats = list(numeric_pat, sne_pat),
                            excl_preds = list(is_all_integer)) {
    frame %>%
        ## filter out column names by regular expression
        Reduce(init = ., x = excl_pats, f = function (df, pat) {
            df %>% select_if({ !grepl(pat, colnames(.), perl = T) })
        }) %>%
        ## filter out columns by their value or type
        Reduce(init = ., x = excl_preds, f = function(df, pred) {
            ## pred takes a vector and outputs a single logical value
            ## apply pred to each column of the data frame and remove columns
            ## for which it returns TRUE
            df %>% select_if({ !apply(., 2, pred) })
        })
}

shared_markers <- function (frames) {
    ## TODO: find less common columns and check if they're mistakes
    ## TODO: if columns are close but not the same (e.g. levenshtein), show a
    ## warning
    frames %>% lapply(colnames) %>% Reduce(f = intersect)
}


### Analyze hierarchies of populations in a dataset.

emd_compute_row <- function (i, n, mats, output, clust, num_cores,
                             max_iterations, verbose) {
    stopifnot(is.integer(n) && (length(n) == 1) &&
              (i <= n) && (length(mats) == n) && all(dim(output) == c(n, n)))
    ## populate result vector and return
    result <- vector('double', n)
    i_mat <- mats[[i]]
    i_rows <- dim(i_mat)[1]
    i_w <- rep(1, i_rows)
    if (i > 1) {
        for (j in 1:(i - 1)) {
            result[j] <- output[j, i]
        }
    }
    result[i] <- 0
    if (i < n) {
        col_range <- (i + 1):n
        if (verbose) {
            cat(sprintf("computing columns %s:%s on %s workers\n",
                        i + 1, n, num_cores))
        }
        ## TODO: see ?parLapply and see if there's an alternative approach
        ## which would be even faster
        parallel::parLapply(clust, col_range, function (j) {
            j_mat <- mats[[j]]
            j_rows <- dim(j_mat)[1]
            j_w <- rep(1, j_rows)
            col_time <- system.time(
                cur_emd <- emdist::emdw(
                    i_mat, i_w,
                    j_mat, j_w,
                    max.iter = max_iterations))
            list(j = j,
                 time = col_time[3],
                 result = cur_emd,
                 worker = Sys.getpid())
        }) %>% lapply(function (r) {
            if (verbose) {
                cat(sprintf("column %s took %s sec on worker %s\n",
                            r$j, r$time, r$worker))
            }
            r$result
        }) %>% unlist -> result[col_range]
    }
    result
}

#' @title Compute pairwise EMD of a set of datasets in a viSNE analysis.
#'
#' @description \code{pairwise_emd} computes the Earth Mover's Distance (EMD)
#'     between the viSNE axes of each pair of data frames given, and writes the
#'     resultant double-precision matrix to a CSV file.
#'
#' @param frames named list of data frames representing CyToF
#'     datasets. \code{\link{process_cyto_dataset}} can be used to produce this.
#' @param outfile string naming a file path or an open connection which the
#'     resultant matrix of comparisons is written to with
#'     \code{\link{write.table}}.
#' @param max_iterations Positive integer indicating the maximum iterations used
#'     to compute EMD. Changing this argument \emph{typically} does not change
#'     the result very much.
#' @param visne_axes character vector indicating the column names containing the
#'     viSNE channels.
#' @param verbose logical. Whether to print progress indicators to the console.
#'
#' @return The value of \code{outfile}. \code{outfile} contains both row and
#'     column names, and should be read back into a matrix for analysis like:
#'     \code{as.matrix(read.csv(outfile, row.names = 1))}.
#'
#' @details If \code{length(frames) == n} for some positive integer \code{n}, a
#'     diagonal n x n matrix \code{mat} is created to represent the EMD between
#'     each pair of datasets at indices \code{i} and \code{j}, where
#'     \code{mat[i,j] == mat[j,i]} and \code{mat[i,j]} represents the EMD
#'     between the two.
#'
#' @seealso \code{\link{emdist::emdw}} for the underlying EMD implementation.
#'
#' @export
#'
pairwise_emd <- function (frames, outfile,
                          max_iterations = 10,
                          visne_axes = c("tSNE1", "tSNE2"),
                          verbose = T) {
    stopifnot(is.vector(visne_axes, 'character') && (length(visne_axes) > 0))
    nm <- get_names(frames)
    n <- length(frames)
    mats <- lapply(frames, function (frame) {
        frame[,visne_axes] %>% as.matrix
    })
    output <- matrix(vector(mode = "double", length = n * n), nrow = n)
    if (verbose) {
        cat(sprintf("starting pairwise emd on %s data frames...\n", n))
    }
    num_cores <- parallel::detectCores()
    clust <- parallel::makeCluster(num_cores)
    tryCatch(finally = parallel::stopCluster(clust), {
        for (i in 1:n) {
            if (verbose) {
                cat(sprintf("row %s/%s\n", i, n))
            }
            row_time <- system.time(
                output[i,] <- emd_compute_row(
                    i, n, mats, output, clust, num_cores, max_iterations,
                    verbose))
            if (verbose) {
                cat(sprintf("row %s/%s took %s seconds to compute %s columns\n",
                            i, n, round(row_time[3], 3), (n - i)))
            }
        }
    })
    colnames(output) <- nm
    rownames(output) <- nm
    write.table(output, outfile, sep = ",")
    stopifnot(file.exists(outfile))
    outfile
}

calc_mag_iqr <- function (frame, markers) {
    lapply(markers, function (mark) {
        frame[,mark] %>% { c(median(.), IQR(., type = 2)) }
    }) %>% Reduce(f = rbind) %>%
        set_colnames(c("MAG", "IQR")) %>%
        set_rownames(markers) %>% as.data.frame
}

calc_mem <- function (pop, ref, markers) {
    stopifnot(all(rownames(pop) == markers) &&
              all(rownames(ref) == markers))
    flip_mems <- (pop$MAG - ref$MAG) < 0
    mems <- (abs(pop$MAG - ref$MAG) + (ref$IQR / pop$IQR) - 1) *
        ((-2 * flip_mems) + 1)
    set_names(mems, markers)
}

#' @title Compute pairwise MEM RMSD of a set of CyToF datasets.
#'
#' @description \code{pairwise_mem_rmsd} computes the Root Mean Square Distance
#'     (RMSD) of the Marker Enrichment Modeling (MEM) score between each each
#'     pair of data frames given, and writes the resultant double-precision
#'     matrix to a CSV file.
#'
#' @param frames named list of data frames representing CyToF
#'     datasets. \code{\link{process_cyto_dataset}} can be used to produce this.
#' @param outfile string naming a file path or an open connection which the
#'     resultant matrix of comparisons is written to with
#'     \code{\link{write.table}}.
#' @param markers character vector of channel names to use for the MEM
#'     calculation, or \code{NULL}. If \code{markers = NULL}, the function will
#'     guess the markers to perform the MEM analysis on. If \code{verbose =
#'     TRUE}, the choice of markers will be printed to the console. If
#'     \code{markers} is non-\code{NULL}, all of the columns indicated should
#'     exist in each input dataset.
#' @param ref_pop data frame containing the reference population to use for the
#'     MEM calculation, or \code{NULL}. If \code{ref_pop = NULL}, each input
#'     dataset will be collapsed (using \code{\link{rbind}}) into a single giant
#'     reference population.
#' @param transform_with specifies how the raw channel values should be
##   transformed for direct comparison.
#' @param verbose logical. Whether to print progress indicators to the console.
#'
#' @return The value of \code{outfile}. \code{outfile} contains both row and
#'     column names, and should be read back into a matrix for analysis like:
#'     \code{as.matrix(read.csv(outfile, row.names = 1))}.
#'
#' @details If \code{length(frames) == n} for some positive integer \code{n}, a
#'     diagonal n x n matrix \code{mat} is created to represent the MEM RMSD
#'     between each pair of datasets at indices \code{i} and \code{j}, where
#'     \code{mat[i,j] == mat[j,i]} and \code{mat[i,j]} represents the MEM RMSD
#'     between the two.
#'
#' @seealso \code{\link{asinh_transform}} is the default transform used to
#'     compare markers.
#'
#' @references Diggins, Kirsten E., et al. "Characterizing cell subsets using
#'     marker enrichment modeling." Nature Methods 14.3 (2017): 275-278.
#'
#' @export
#'
pairwise_mem_rmsd <- function (frames, outfile,
                          markers = NULL, ref_pop = NULL,
                          transform_with = asinh_transform,
                          verbose = T) {
    nm <- get_names(frames)
    n <- length(frames)
    if (verbose) {
        cat(sprintf("performing pairwise MEM on %s data frames...\n", n))
    }
    on_markers <-
        if (!is.null(markers)) {
            markers
        } else {
            if (verbose) {
                cat("guessing markers to join on...\n")
            }
            frames %>% lapply(cyto_data_cols) %>% shared_markers
        }
    if (verbose) {
        cat(sprintf("joining on markers:\n[%s]\n",
                    paste0(on_markers, collapse = ", ")))
        cat("performing MEM...\n")
    }
    marked_pops <- lapply(frames, function (df) {
        df[,on_markers] %>% mutate_all(transform_with)
    })
    stat_dfs <- lapply(marked_pops, function (df) calc_mag_iqr(df, on_markers))
    used_ref_pop <- if (!is.null(ref_pop)) {
                        ref_pop
                    } else { Reduce(x = marked_pops, f = rbind) }
    ref_mem <- calc_mag_iqr(used_ref_pop, on_markers)
    mem_vectors <- lapply(stat_dfs, function (st_df) {
        calc_mem(st_df, ref_mem, on_markers)
    })
    output <- matrix(vector("double", length = n * n), nrow = n)
    if (verbose) {
        cat(sprintf("starting pairwise MEM RMSD on %s data frames...\n", n))
    }
    for (i in 1:n) {
        if (verbose) {
            cat(sprintf("row %s/%s\n", i, n))
        }
        i_vt <- mem_vectors[[i]]
        for (j in 1:n) {
            if (verbose) {
                cat(sprintf("col %s/%s\n", j, n))
            }
            j_vt <- mem_vectors[[j]]
            output[i,j] <- (i_vt - j_vt) ^ 2 %>% sum %>% sqrt
        }
    }
    colnames(output) <- nm
    rownames(output) <- nm
    write.table(output, outfile, sep = ",")
    stopifnot(file.exists(outfile))
    outfile
}

#' @title Plot Pairwise Comparisons of Cytometry Datasets
#'
#' @description \code{plot_pairwise_comparison} plots a heatmap of a pairwise
#'     comparison of datasets with \code{\link{gplots::heatmap.2}}.
#'
#' @param matrix_file string naming a file path or an open connection which a
#'     comparison matrix is read from. File should be produced by
#'     \code{\link{pairwise_emd}} or \code{\link{pairwise_mem_rmsd}}.
#' @param color_palette color palette generated by
#'     \code{\link{grDevices::colorRampPalette}}, or \code{NULL}.
#'
#' @seealso \code{\link{gplots::heatmap.2}} for the underlying plotting
#'     function, and \code{\link{grDevices::colorRampPalette}} for color palette
#'     generation. \code{\link{pairwise_emd}} or \code{\link{pairwise_mem_rmsd}}
#'     should be used to generate \code{matrix_file}.
#'
#' @export
#'
plot_pairwise_comparison <- function (matrix_file, color_palette = NULL) {
    mat <- as.matrix(read.csv(matrix_file, row.names = 1))
    if (is.null(color_palette)) {
        gplots::heatmap.2(
            mat,
            Rowv = F, Colv = F, dendrogram = "none",
            trace = "none", density.info = "none")
    } else {
        gplots::heatmap.2(
            mat, col = color_palette,
            Rowv = F, Colv = F, dendrogram = "none",
            trace = "none", density.info = "none")
    }
}


### Wrapper functions for xpath selectors.

xpath <- function (doc, xpath_str = ".", node = doc, fun = NULL) {
    doc_ns <- XML::xmlNamespaceDefinitions(doc, simplify = T)
    XML::xpathSApply(
        doc = node, path = xpath_str, fun = fun, namespaces = doc_ns)
}



### Parsing cytobank gatingML xml files.

## gen_generic <- function (name, ) {
## }

setGeneric("add_gate", function (gate, cyto_df) {

})

setClass("Gate",
         slots = c(gate_name = "character", gate_id = "character"))
setMethod("initialize", "Gate", function(.Object, gate_xml, doc) {
    .Object <- callNextMethod(.Object)
    gate_name <- xpath(doc, "data-type:custom_info/cytobank/name/text()",
                       gate_xml, xmlValue) %>% get_single
    gate_id <- xmlGetAttr(gate_xml, "gating:id")
    .Object@gate_name <- gate_name
    .Object@gate_id <- gate_id
    .Object
})

data_transformation_fun_dict <- list(
    Tr_Arcsinh_5 = asinh_transform
)

setClass("RectangleGate",
         slots = c(constraints = "list"),
         contains = "Gate")
setMethod(
    "initialize", "RectangleGate",
    function (.Object, rect_gate_xml, doc) {
        .Object <- callNextMethod(.Object, rect_gate_xml, doc)
        constraints <- xpath(doc, "gating:dimension", rect_gate_xml) %>%
            lapply(function (dim_xml) { list(
                transform_name = xpath(
                    doc, "@gating:transformation-ref", dim_xml),
                marker = xpath(
                    doc, "data-type:fcs-dimension/@data-type:name", dim_xml),
                min = xpath(doc, "@gating:min", dim_xml),
                max = xpath(doc, "gating:max", dim_xml))
            })
        .Object@constraints <- constraints
        .Object
    })
setMethod(
    "add_gate", c(gate = "RectangleGate", cyto_df = "data.frame"),
    function (gate, cyto_df) {
        cyto_df[,gate@gate_name] <- TRUE
        Reduce(
            init = cyto_df,
            x = gate@constraints,
            f = function (acc, constr) {
                tr_f <- match.fun(
                    data_transformation_fun_dict[[constr$transform_name]])
                satisfies_constraint <- acc[,constr$marker] %>% tr_f %>% {
                    (. >= constr$min) & (. <= constr$max)
                }
                acc[,gate@gate_name] <-
                    acc[,gate@gate_name] & satisfies_constraint
                acc
            })
    })

## example:
## https://github.com/RGLab/CytoML/blob/trunk/R/read.gatingML.cytobank.R

## setClassUnion(
##     "Gate", c("RectangleGate", "PolygonGate", "BooleanGate", "QuadrantGate"))

## doc <- xmlParse("./allie-paper/CytExp_22899_Gates_v1.xml")
## cyto_df <- list.files(
##     path = "./allie-paper/", pattern = "fcs$", full.names = T)
## gates <- parse_rectangle_gates(doc)
## processed_cyto_frames <- lapply(cyto_df, function (cyto_data_file) {
##     Reduce(init = read_file(cyto_data_file), x = gates,
##            f = function (acc, cur) {
##                add_gate(cur, acc)
##            })
## })

parse_rectangle_gates <- function (xml) {
    ## "data-type:custom_info/cytobank/fcs_file_filename/text()"
    xpath(xml, "/gating:Gating-ML/gating:RectangleGate") %>%
        lapply(function (rect_gate_node) {
            new("RectangleGate", rect_gate_node, xml)
        })
}

parse_polygon_gates <- function (xml) {
    xml
}

parse_quadrant_gates <- function (xml) {
    quadrant_gate_nodes <- xpath(xml, "/gating:Gating-ML/gating:QuadrantGate")
    lapply(quadrant_gate_nodes, function (quadrant_gate_node) {
    })
}

parse_boolean_gates <- function (xml) {
    xml
}

gate_parse_dispatch <- list(
    ## TODO: check if there are node names matching /^gating:.*Gate$/ that
    ## aren't in this list
    PolygonGate = parse_polygon_gates,
    RectangleGate = parse_rectangle_gates,
    QuadrantGate = parse_quadrant_gates,
    BooleanGate = parse_boolean_gates
)
