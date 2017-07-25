#' @import magrittr
#' @import dplyr


### Clean fcs data.

#' @description ?
#'
#' @export
#' @rdname non_pheno_channel_name_patterns
#'
numeric_channel_name <- "^[[:digit:]]+$"
#' @description ?
#'
#' @export
#' @rdname non_pheno_channel_name_patterns
#'
sne_channel_name <- "[sS][nN][eE]"
#' @title ?
#'
#' @description ?
#'
#' @export
#'
non_pheno_channel_name_patterns <- c(
    numeric_channel_name,
    sne_channel_name)

## FIXME: does read.FCS read in anything except double columns?
#' @description ?
#'
#' @param vec ?
#'
#' @return ?
#'
#' @export
#' @rdname non_pheno_channel_content_predicates
#'
is_non_floating_point <- function (vec) {
    (!is.vector(vec, mode = 'double')) ||
        all(vec == as.integer(vec))
}
## FIXME: check for marker columns with "spread" of data as heuristic
## (e.g. viSNE columns are between +/-30)
#' @title ?
#'
#' @description ?
#'
#' @export
#'
non_pheno_channel_content_predicates <- c(is_non_floating_point)

#' @title ?
#'
#' @description ?
#'
#' @param frames ?
#' @param name_filter ?
#' @param content_filter ?
#' @param check.spelling ?
#' @param check.channel_switches ?
#'
#' @return ?
#'
#' @details ?
#'
#' @export
#'
get_marker_names <- function
(
    frames,
    name_filter = non_pheno_channel_name_patterns,
    content_filter = non_pheno_channel_content_predicates,
    check.spelling = TRUE,
    check.channel_switches = TRUE
) {
    nm <- get_names(frames)
    frames_filtered_cols <- lapply(frames, function (fcs_df) {
        colnames(fcs_df) %>%
            Reduce(init = ., x = name_filter, f = function (cols, pat) {
                removed_colnames <- grepl(pat, cols, perl = TRUE)
                cols[!removed_colnames]
            }) %>%
            select(fcs_df, .)
    })
    if (check.spelling) {
        cols_by_frame <- lapply(frames_filtered_cols, colnames)
        all_cols <- Reduce(f = c, x = cols_by_frame) %>% unique
        lapply(cols_by_frame, function (cur_cols) {

        })
        lapply(frames_filtered_cols, function (fcs_df) {

        })
    }

    ## lapply(frames, function (fcs_df) {
    ##     ## filter out column names by regular expressions
    ##     Reduce(
    ##         init = colnames(fcs_df), x = name_filter,
    ##         f = function (cols, pat) {
    ##             removed_colnames <- grepl(pat, cols, perl = TRUE)
    ##             cols[!removed_colnames]
    ##         })
    ##     fcs_df[,]
    ## })
    ## filtered_frames <- lapply(frames, function (fcs_df) {
    ##     fcs_filtered_colnames <- fcs_df[,removed_colnames]
    ##     ## filter out columns by functions of their content
    ##     fcs_markers_only <- Reduce(
    ##         init = fcs_filtered_colnames, x = content_filter,
    ##         f = function (remaining_df, pred) {
    ##             removed_cols <- apply(remaining_df, 2, pred)
    ##             remaining_df[,!removed_cols]
    ##         })
    ##     fcs_markers_only
    ## })
    ## if (check.spelling) {
    ## }
}

## Return a char vector containing the column names shared among all data frames
## given.
shared_markers <- function (frames) {
    ## FIXME: find less common columns and check if they're mistakes
    ## FIXME: if columns are close but not the same (e.g. levenshtein), show a
    ## warning
    ## TODO: offer option for MEM RMSD on the marker names shared between each
    ## pair of files, not just the columns shared between ALL files (and ENSURE
    ## this is noted and explained)
    frames %>% lapply(colnames) %>% Reduce(f = intersect)
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
#' @param frames named list of data frames representing CyToF datasets.
#' @param visne_axes character vector indicating the column names containing the
#'     viSNE channels.
#' @param verbose logical. Whether to print progress indicators to the console.
#'
#' @return The resultant matrix of EMD between all pairs of input files's viSNE
#'     axis values. Row and column names are set to the input file paths.
#'
#' @details If \code{length(frames) == n} for some positive integer \code{n}, a
#'     diagonal n x n double-precision floating-point matrix \code{mat} is
#'     created to represent the EMD between each pair of datasets at indices
#'     \code{i} and \code{j}, where \code{mat[i,j] == mat[j,i]} and
#'     \code{mat[i,j]} represents the EMD between the two.
#'
#' @seealso \code{\link{transport::wasserstein}} for the underlying EMD
#'     implementation.
#'
#' @export
#'
pairwise_emd <- function (frames,
                          visne_axes = c("tSNE1", "tSNE2"),
                          verbose = TRUE) {
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
    output %>% set_colnames(nm) %>% set_rownames(nm)
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
#'     datasets.
#' @param markers character vector of channel names to use for the MEM
#'     calculation. All of the columns indicated should exist in each input
#'     dataset.
#' @param ref_pop data frame produced by \code{\link{read_cyto_file}} containing
#'     the reference population to use for the MEM calculation, or
#'     \code{NULL}. If \code{ref_pop = NULL}, each input dataset will be
#'     collapsed (using \code{\link{rbind}}) into a single giant reference
#'     population.
#' @param transform_with specifies how the raw channel values should be
#'     transformed for direct comparison.
#' @param verbose logical. Whether to print progress indicators to the console.
#'
#' @return The resultant matrix of MEM RMSD between all pairs of input files
#'     along the specified \code{markers}. Row and column names are set to the
#'     input file paths.
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
pairwise_mem_rmsd <- function (frames, markers,
                               ref_pop = NULL,
                               transform_with = asinh_transform,
                               verbose = TRUE) {
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
    output %>% set_colnames(nm) %>% set_rownames(nm)
}


#' @title Plot Pairwise Comparisons of Cytometry Datasets
#'
#' @description \code{plot_pairwise_comparison} plots a heatmap of a pairwise
#'     comparison of datasets with \code{\link{gplots::heatmap.2}}.
#'
#' @param mat A comparison matrix produced by \code{\link{pairwise_emd}} or
#'     \code{\link{pairwise_mem_rmsd}}.
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
plot_pairwise_comparison <- function (mat,
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
