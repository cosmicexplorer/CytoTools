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
non_pheno_channel_name_patterns <- list(
    "^[[:digit:]]+$" = NA,
    "[sS][nN][eE]" = NA,
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
    frame, desc,
    content_filter = non_pheno_channel_content_predicates,
    normalize_channel_names = toupper,
    channel_name_ops = non_pheno_channel_name_patterns
) {
    by_content <- Reduce(
        init = frame, x = content_filter,
        f = function (df, pred) select_if(df, Negate(match.fun(pred))))
    colnames(by_content) %>%
        { match.fun(normalize_channel_names)(.) } %>%
        replace_matches(channel_name_ops) %>%
        replace_colnames(by_content, desc, .)
}

## TODO: try to find non-marker columns with "spread" of data as heuristic
## (e.g. viSNE columns are between +/-30) and throw if any found
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
        ## NOTE: assignment with -> at bottom!
        parallel::parLapply(clust, col_range, function (j) {
            j_wpp <- measures[[j]]
            ## TODO: allow controlling parameters of this calculation!
            col_time <- system.time(
                gcFirst = FALSE,
                cur_emd <- transport::wasserstein(
                    i_wpp, j_wpp, control = transport::trcontrol(
                        a = i_wpp, b = j_wpp, nscales = 3, scmult = 3)))
            ## can't use msg() here because within worker
            message(sprintf("column %s took %s sec on worker id %s",
                            j, round(col_time[3], 3), Sys.getpid()))
            cur_emd
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
    msg("starting pairwise emd on %s populations...", n)
    num_cores <- parallel::detectCores()
    ## NOTE: outfile = "" goes to console, EXCEPT in Rgui on Windows!
    clust <- parallel::makeCluster(num_cores, outfile = "")
    time_tot <- system.time(
        gcFirst = FALSE,
        tryCatch(finally = parallel::stopCluster(clust), {
            for (i in 1:n) {
                msg("row %s/%s", i, n)
                row_time <- system.time(
                    gcFirst = FALSE,
                    output[i,] <- emd_compute_row(
                        i, n, measures, output, clust, num_cores))
                msg("row %s/%s took %s seconds to compute %s columns",
                    i, n, round(row_time[3], 3), (n - i))
            }
        }))
    msg("pairwise emd on %s populations took %s seconds",
        n, round(time_tot[3], 3))
    output
}

calc_mag_iqr <- function (pop) {
    apply(pop, 2, function (col) {
        c(median(col), IQR(col, type = 2))
    }) %>% set_rownames(c("MAG", "IQR")) %>% t %>% as.data.frame
}

calc_mem <- function (pop, ref) {
    pop_markers <- rownames(pop)
    rownames(ref) %>% {
        stopifnot(!is.null(pop_markers) && !is.null(.) &&
                  all(pop_markers == .))
    }
    flip_factors <- (-2 * (pop$MAG < ref$MAG)) + 1
    mems <- flip_factors * (abs(pop$MAG - ref$MAG) + (ref$IQR / pop$IQR) - 1)
    setNames(mems, pop_markers)
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
#' @param squash_fun specifies how the raw channel values should be
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
pairwise_mem_rmsd <- function (pops,
                               ref_pop = NULL,
                               squash_fun = asinh_transform,
                               ...) {
    nm <- get_names(pops)
    n <- length(pops)
    normalized_pops <- lapply(1:n, function (i) {
        normalize_channels(pops[[i]], nm[[i]], ...)
    })
    shared_channels <- normalized_pops %>%
        lapply(colnames) %>% Reduce(f = intersect)
    norm_ref_pop <- if (!is.null(ref_pop)) {
                        normalize_channels(ref_pop, ...)
                    } else {
                        ## don't need to norm again
                        normalized_pops %>%
                            lapply((. %>% select(shared_channels))) %>%
                            Reduce(f = rbind)
                    }
    ## noop if ref_pop is null
    shared_channels <- intersect(shared_channels, colnames(norm_ref_pop))
    msg("joining on %s markers: [%s]",
        length(shared_channels), paste0(shared_channels, collapse = ", "))
    msg("performing pairwise MEM on %s populations...", n)
    ## TODO: check channels here and throw if misspelling/etc
    get_pop_stats <-
        (. %>%
         select(shared_channels) %>%
         mutate_all(match.fun(squash_fun)) %>%
         calc_mag_iqr)
    stats_by_pop <- normalized_pops %>% lapply(get_pop_stats)
    ref_stats <- get_pop_stats(norm_ref_pop)
    mems_by_pop <- stats_by_pop %>% lapply((. %>% calc_mem(., ref_stats)))
    output <- matrix(double(n * n), n, n, dimnames = list(nm, nm))
    msg("starting pairwise MEM RMSD on %s data frames...", n)
    for (i in 1:n) {
        msg("row %s/%s", i, n)
        i_mem <- mems_by_pop[[i]]
        for (j in 1:n) {
            msg("col %s/%s", j, n)
            j_mem <- mems_by_pop[[j]]
            output[i,j] <- ((i_mem - j_mem) ^ 2) %>% sum %>% sqrt
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
plot_pairwise_comparison <- function (mat, dendro = FALSE, ...) {
    opts <- merge_named_lists(
        list(...),
        list(trace = "none", density.info = "none",
             cexRow = .5, cexCol = .5, margin = c(10, 10)))
    if (!dendro) {
        opts <- merge_named_lists(
            opts,
            list(Rowv = FALSE,
                 Colv = FALSE,
                 dendrogram = "none"))
    }
    map_fun <- get("heatmap.2", asNamespace("gplots"))
    arglist <- c(list(mat), opts)
    do.call(map_fun, arglist)
}
