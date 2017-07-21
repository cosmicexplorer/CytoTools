#' @import magrittr
#' @import dplyr


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

## Return data frame which contains (probably) only marker data.
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

## Return a char vector containing the column names shared among all data frames
## given.
shared_markers <- function (frames) {
    ## TODO: find less common columns and check if they're mistakes
    ## TODO: if columns are close but not the same (e.g. levenshtein), show a
    ## warning
    frames %>% lapply(colnames) %>% Reduce(f = intersect)
}


### Analyze hierarchies of populations in a dataset.

emd_compute_row <- function (i, n, measures, output, clust, num_cores,
                             verbose,
                             nscales = 3, scmult = 3, ...) {
    ## only used for argument checking
    cmb <- c(i, n)
    stopifnot(is.integer(cmb) &&
              (length(cmb) == 2) &&
              (i <= n) &&
              (length(measures) == n) &&
              all(dim(output) == c(n, n)))
    ## populate result vector and return
    result <- vector('double', n)
    i_wpp <- measures[[i]]
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
            j_wpp <- measures[[j]]
            ## TODO: allow controlling iterations somehow!
            col_time <- system.time(
                cur_emd <- transport::wasserstein(
                    i_wpp, j_wpp, control = trcontrol(
                        a = i_wpp, b = j_wpp,
                        nscales = nscales, scmult = scmult, ...)))
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
                          visne_axes = c("tSNE1", "tSNE2"),
                          verbose = T) {
    stopifnot(is.vector(visne_axes, 'character') && (length(visne_axes) > 0))
    nm <- get_names(frames)
    n <- length(frames)
    measures <- lapply(frames, function (frame) {
        mat <- frame[,visne_axes] %>% as.matrix
        transport::wpp(mat, rep(1, dim(mat)[1]))
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
                    i, n, measures, output, clust, num_cores, verbose))
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
#' @param ref_pop data frame produced by \code{\link{read_cyto_file}} containing
#'     the reference population to use for the MEM calculation, or
#'     \code{NULL}. If \code{ref_pop = NULL}, each input dataset will be
#'     collapsed (using \code{\link{rbind}}) into a single giant reference
#'     population.
#' @param transform_with specifies how the raw channel values should be #
#transformed for direct comparison.
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
#' @param dendro logical. whether to display a dendrogram in the output plot.
#'
#' @seealso \code{\link{gplots::heatmap.2}} for the underlying plotting
#'     function, and \code{\link{grDevices::colorRampPalette}} for color palette
#'     generation. \code{\link{pairwise_emd}} or \code{\link{pairwise_mem_rmsd}}
#'     should be used to generate \code{matrix_file}.
#'
#' @export
#'
plot_pairwise_comparison <- function (matrix_file,
                                      color_palette = NULL,
                                      dendro = FALSE) {
    mat <- as.matrix(read.csv(matrix_file, row.names = 1))
    arglist <- list(mat, trace = "none", density.info = "none")
    if (!is.null(color_palette)) {
        arglist <- c(arglist, list(col = color_palette))
    }
    if (!dendro) {
        arglist <- c(arglist, list(Rowv = FALSE,
                                   Colv = FALSE,
                                   dendrogram = "none"))
    }
    map_fun <- get("heatmap.2", asNamespace("gplots"))
    eval(do.call(map_fun, arglist))
}
