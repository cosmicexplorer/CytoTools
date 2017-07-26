#' @import magrittr
#' @import dplyr


### Clean fcs data.

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
#' @export
#'
non_pheno_channel_name_patterns <- c(
    "^[[:digit:]]+$" = NA_character_,
    "TSNE" = NA_character_,
    "\\-[0-9]+$" = "")

#' @title ?
#'
#' @description ?
#'
#' @param frame ?
#'
#' @return ?
#'
#' @export
#'
normalize_channels <- function
(
    frame,
    content_filter = non_pheno_channel_content_predicates,
    normalize_channel_names = toupper,
    channel_name_ops = non_pheno_channel_name_patterns
) {
    by_content <- Reduce(
        init = frame, x = content_filter,
        f = function (df, pred) select_if(df, match.fun(pred)))
    colnames(by_content) %>%
        { match.fun(normalize_channel_names)(.) } %>%
        str_replace_all(channel_name_ops) %>%
        replace_colnames(by_content, .)
}

check_marker_spelling <- function (frames
                                 ## , maxdist, norm = tolower
                                   ) {
    ## cols_by_frame <- lapply(frames, colnames)
    ## shared_cols <- Reduce(f = intersect, x = cols_by_frame)
    ## lapply(cols_by_frame, function (cur_cols) {
    ##     cur_cols
    ## })
}

check_channel_switches <- function (frames) {
}

## TODO: try to find non-marker columns with "spread" of data as heuristic
## (e.g. viSNE columns are between +/-30) and throw if any found
#' @title ?
#'
#' @description ?
#'
#' @param frames ?
#'
#' @return ?
#'
#' @details ?
#'
#' @export
#'
normalize_channels <- function (frames) {
    nm <- get_names(frames)
    frames_filtered_cols
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

emd_compute_row <- function (i, n, measures, output, clust, num_cores) {
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
        msg("computing columns %s:%s on %s workers", i + 1, n, num_cores)
        ## TODO: see ?parLapply and see if there's an alternative approach
        ## which would be even faster
        parallel::parLapply(clust, col_range, function (j) {
            j_wpp <- measures[[j]]
            ## TODO: allow controlling parameters of this calculation!
            col_time <- system.time(
                cur_emd <- transport::wasserstein(
                    i_wpp, j_wpp, control = transport::trcontrol(
                        a = i_wpp, b = j_wpp, nscales = 3, scmult = 3)))
            list(j = j,
                 time = col_time[3],
                 result = cur_emd,
                 worker = Sys.getpid())
        }) %>% lapply(function (r) {
            msg("column %s took %s sec on worker %s", r$j, r$time, r$worker)
            r$result
        }) %>% unlist -> result[col_range]
    }
    result
}


#' @title Compute pairwise EMD of a set of datasets in a viSNE analysis.
#'
#' @description ? \code{pairwise_emd} computes the Earth Mover's Distance (EMD)
#'     between the viSNE axes of each pair of data frames given, and writes the
#'     resultant double-precision matrix to a CSV file.
#'
#' @param tsne_matrices ?
#'
#' @return ? The resultant matrix of EMD between all pairs of input files's
#'     viSNE axis values. Row and column names are set to the input file paths.
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
pairwise_emd <- function (tsne_matrices) {
    nm <- get_names(tsne_matrices)
    n <- length(tsne_matrices)
    measures <- lapply(tsne_matrices, function (mat) {
        transport::wpp(mat, rep(1, dim(mat)[1]))
    })
    output <- matrix(double(n * n), n, n, dimnames = list(nm, nm))
    msg("starting pairwise emd on %s data frames...", n)
    num_cores <- parallel::detectCores()
    clust <- parallel::makeCluster(num_cores)
    tryCatch(finally = parallel::stopCluster(clust), {
        for (i in 1:n) {
            msg("row %s/%s", i, n)
            row_time <- system.time(
                output[i,] <- emd_compute_row(
                    i, n, measures, output, clust, num_cores))
            msg("row %s/%s took %s seconds to compute %s columns",
                i, n, round(row_time[3], 3), (n - i))
        }
        output
    })
}

calc_mag_iqr <- function (frame) {
    apply(frame, 2, function (col) {
        c(median(col), IQR(col, type = 2))
    }) %>% set_rownames(c("MAG", "IQR")) %>% t
}

calc_mem <- function (pop, ref) {
    flip_mems <- pop$MAG < ref$MAG
    flip_factors <- (-2 * flip_mems) + 1
    flip_factors * (abs(pop$MAG - ref$MAG) + (ref$IQR / pop$IQR) - 1)
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
pairwise_mem_rmsd <- function (frames,
                               ref_pop = NULL,
                               transform_with = asinh_transform) {
    stopifnot(is.list(frames) && length(frames) >= 1)
    
    markers <- colnames(frames[[1]])
    lapply(frames, function (fcs_df) {
        stopifnot(setequal(colnames(fcs_df), markers))
    })
    nm <- get_names(frames)
    n <- length(frames)
    msg("performing pairwise MEM on %s data frames...", n)
    marked_pops <- lapply(frames, function (df) {
        df[,markers] %>% mutate_all(match.fun(transform_with))
    })
    stat_dfs <- lapply(marked_pops, calc_mag_iqr)
    used_ref_pop <- if (!is.null(ref_pop)) {
                        ref_pop
                    } else { Reduce(x = marked_pops, f = rbind) }
    ref_mem <- calc_mag_iqr(used_ref_pop)
    mem_vectors <- lapply(stat_dfs, function (st_df) {
        calc_mem(st_df, ref_mem) %>% set_names(markers)
    })
    output <- matrix(double(n * n), n, n, dimnames = c(nm, nm))
    msg("starting pairwise MEM RMSD on %s data frames...", n)
    for (i in 1:n) {
        msg("row %s/%s", i, n)
        i_vt <- mem_vectors[[i]]
        for (j in 1:n) {
            msg("col %s/%s", j, n)
            j_vt <- mem_vectors[[j]]
            output[i,j] <- (i_vt - j_vt) ^ 2 %>% sum %>% sqrt
        }
    }
    output
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
