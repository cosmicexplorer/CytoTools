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
#' @param pop_list ?
#' @param ref ?
#' @param transform_fun ?
#' @param content_filter ?
#' @param modify_colnames_fun ?
#' @param channel_name_ops ?
#'
#' @return ?
#'
#' @export
#'
normalize_pheno_channels_dataset <- function
(
    pop_list,
    ref = NULL,
    transform_fun = asinh_transform,
    content_filter = non_pheno_channel_content_predicates,
    modify_colnames_fun = toupper,
    channel_name_ops = non_pheno_channel_name_patterns
) {
    gen_modify_fun <- function (pop, desc) {
        by_content <- Reduce(
            init = pop, x = content_filter,
            f = function (df, pred) select_if(df, Negate(pred)))
        colnames(by_content) %>%
            modify_colnames_fun %>%
            replace_matches(channel_name_ops) %>%
            replace_colnames(by_content, desc, .) %>%
            as.matrix %>%
            transform_fun
    }
    stopifnot(is.list(pop_list))
    pop_names <- names(pop_list)
    indiv_normed_channel_pops <- lapply(1:length(pop_list), function (i) {
        gen_modify_fun(pop_list[[i]], pop_names[[i]])
    }) %>% set_names(pop_names)
    shared_channels <- indiv_normed_channel_pops %>%
        lapply(colnames) %>% Reduce(f = intersect)
    if (!is.null(ref)) {
        norm_ref_pop <- gen_modify_fun(ref, "REF")
        shared_channels <- intersect(shared_channels, colnames(norm_ref_pop))
    } else {
        norm_ref_pop <- indiv_normed_channel_pops %>%
            lapply((. %>% .[,shared_channels])) %>%
            Reduce(rbind)
    }
    list(shared_channels = shared_channels,
         ref = norm_ref_pop[,shared_channels],
         pop_list = indiv_normed_channel_pops %>%
             lapply((. %>% .[,shared_channels])))
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


### Analyze populations in a dataset.

compute_pairwise_emd_trials <- function (matrices, clust,
                                         downsample_rows, comparison_runs,
                                         verbose_timing) {
    n <- length(matrices)
    msg("starting pairwise emd on %s populations...", n)
    lapply(1:n, function (i) {
        final_col <- as.integer(n - i)
        msg("computing row %s/%s, requiring column(s) %s:%s",
            i, n, 1, final_col)
        i_mat <- matrices[[i]]
        cols <- safe_int_seq(1L, final_col)
        parallel::parLapply(clust, cols, function (j) {
            this_pid <- Sys.getpid()
            j_mat <- matrices[[n - j + 1]]
            lapply(1:comparison_runs, function (r) {
                i_wpp <- sample(1:nrow(i_mat), downsample_rows) %>%
                    i_mat[.,] %>% transport::wpp(rep(1, downsample_rows))
                j_wpp <- sample(1:nrow(j_mat), downsample_rows) %>%
                    j_mat[.,] %>% transport::wpp(rep(1, downsample_rows))
                single_entry <- CytoTools::timed_execute(transport::wasserstein(
                        i_wpp, j_wpp, control = transport::trcontrol(
                            a = i_wpp, b = j_wpp, method = "shortsimplex")),
                        ## weird that this is necessary when using timed_execute
                        ## inside a worker
                        env = environment())
                if (verbose_timing) {
                    CytoTools::msg(paste(
                        "run %s/%s for row %s/%s, column %s/%s",
                        "on worker %s took %s sec"),
                        r, comparison_runs, i, n, j, length(cols),
                        this_pid, single_entry$time)
                }
                single_entry$value
            }) %>% unlist
        })
    })
}

summarize_computed_trials <- function (fun, trials, nm) {
    stopifnot(is.function(fun))
    n <- length(trials)
    result <- matrix(double(n * n), n, n, dimnames = list(nm, rev(nm)))
    for (i in 1:n) {
        next_cols <- result[rev(safe_int_seq(1L, as.integer(i - 1))),
                            (n - i + 1)]
        result[i,] <- lapply(trials[[i]], fun) %>% unlist %>% c(0, next_cols)
    }
    result
}

#' @title ?
#'
#' @description ?
#'
#' @export
#'
emd_downsample_rows_default <- 1000L
#' @title ?
#'
#' @description ?
#'
#' @export
#' @rdname emd_downsample_rows_default
#'
emd_comparison_runs_default <- 10L


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
pairwise_emd <- function (tsne_matrices,
                          downsample_rows = emd_downsample_rows_default,
                          comparison_runs = emd_comparison_runs_default,
                          summary_funs = list(mean = mean,
                                              variance = var),
                          parallel_cores = parallel::detectCores(),
                          verbose_timing = FALSE) {
    stopifnot(is_counting_num(downsample_rows),
              is_counting_num(comparison_runs),
              is_counting_num(parallel_cores))
    nm <- get_names(tsne_matrices)
    ## NOTE: outfile = "" goes to console, EXCEPT in Rgui on Windows!
    clust <- parallel::makeCluster(parallel_cores, outfile = "")
    computed_trials <- tryCatch(
        finally = parallel::stopCluster(clust),
        compute_pairwise_emd_trials(
            tsne_matrices, clust, downsample_rows, comparison_runs,
            verbose_timing))
    c(list(downsample_rows = downsample_rows,
           comparison_runs = comparison_runs),
      lapply(summary_funs, (. %>% summarize_computed_trials(
          computed_trials, nm))))
}

normal_iqr <- function (x, ...) {
    IQR(x, type = 2, ...)
}

#' @title ?
#'
#' @description ?
#'
#' @references Diggins, Kirsten E., et al. "Characterizing cell subsets using
#'     marker enrichment modeling." Nature Methods 14.3 (2017): 275-278.
#'
#' @export
#'
mem_iqr_threshold_default <- 0.5

#' @title ?
#'
#' @description
#'
#' @param pop_list ?
#' @param ref ?
#' @param IQRthresh ?
#' @param scale_limit ?
#'
#' @return ?
#'
#' @details ?
#'
#' @references Diggins, Kirsten E., et al. "Characterizing cell subsets using
#'     marker enrichment modeling." Nature Methods 14.3 (2017): 275-278.
#'
#' @export
#'
calc_mem <- function (pop_list, ref,
                      IQRthresh = mem_iqr_threshold_default,
                      scale_limit = NULL) {
    markers <- colnames(ref)
    for (pop in pop_list) {
        stopifnot(compare_names(colnames(pop), markers))
    }
    pop_names <- names(pop_list)
    ## make data frames of pop median/iqr
    mag_pops <- pop_list %>%
        lapply((. %>% apply(2, median))) %>%
        Reduce(f = cbind) %>%
        t %>% set_rownames(pop_names)
    iqr_pops <- pop_list %>%
        lapply((. %>% apply(2, normal_iqr))) %>%
        Reduce(f = cbind) %>%
        t %>% set_rownames(pop_names)
    ## vectors for ref stats
    MAGref <- ref %>% apply(2, median) %>% t
    IQRref <- ref %>% apply(2, normal_iqr) %>% t
    ## FIXME: current MEM function takes abs of MAGpop and MAGref (WHY???)
    mag_diffs <- mag_pops %>%
        apply(1, function (MAGpop) abs(MAGpop) - abs(MAGref)) %>%
        t %>% set_colnames(markers)
    ## FIXME: current MEM function thresholds IQRref as well (WHY???)
    stopifnot(is.vector(IQRthresh, mode = 'double'), IQRthresh > 0)
    iqr_ratios_with_threshold <- pmax(iqr_pops, IQRthresh) %>%
        apply(1, function (IQRpop) pmax(IQRref, IQRthresh) / IQRpop) %>%
        t %>% set_colnames(markers)
    ## this is the MEM formula -- see reference
    mems <- (abs(mag_diffs) + iqr_ratios_with_threshold - 1) * sign(mag_diffs)
    ## dilate so the greatest MEM value has magnitude equal to scale_limit
    if (is.null(scale_limit)) {
        mems
    } else {
        max_mem <- mems %>% abs %>% max
        mems / max_mem * scale_limit
    }
}


#' @title Plot Pairwise Comparisons of Cytometry Datasets
#'
#' @description \code{plot_pairwise_comparison} plots a heatmap of a pairwise
#'     comparison of datasets with \code{\link{gplots::heatmap.2}}.
#'
#' @param mat A comparison matrix produced by \code{\link{pairwise_emd}} or
#'     \code{\link{pairwise_mem_rmsd}}.
#' @param is_dist ?
#' @param with_dendrograms ?
#' @param ... ?
#'
#' @seealso \code{\link{gplots::heatmap.2}} for the underlying plotting
#'     function, and \code{\link{grDevices::colorRampPalette}} for color palette
#'     generation. \code{\link{pairwise_emd}} or \code{\link{pairwise_mem_rmsd}}
#'     should be used to generate \code{matrix_file}.
#'
#' @export
#'
plot_pairwise_comparison <- function (mat,
                                      is_dist = TRUE,
                                      no_rotate = FALSE,
                                      with_dendrograms = FALSE,
                                      ...) {
    opts <- merge_named_lists(
        list(...),
        list(trace = "none", scale = "none",
             cexRow = .5, cexCol = .5, margin = c(10, 10)))
    if (is_dist) {
        dims <- dim(mat)
        stopifnot(dims[1] == dims[2],
                  all(mat >= 0))
        ## dist is >= 0 -- this scales it from 0 - 100, where 0 is max dist, 100
        ## is 0 dist
        mat_greatest <- mat %>% max
        mat <- 100 - (mat / mat_greatest * 100)
        opts <- merge_named_lists(opts, list(
            density.info = "histogram",
            symm = TRUE))
    }
    if (!with_dendrograms) {
        opts <- merge_named_lists(opts, list(
            Rowv = FALSE, Colv = FALSE, dendrogram = "none"))
    }
    map_fun <- get("heatmap.2", asNamespace("gplots"))
    arglist <- c(list(mat), opts)
    do.call(map_fun, arglist)
}
