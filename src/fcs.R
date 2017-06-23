### Routines for manipulating CyToF data and interacting with Cytobank.
## Written by Danny McClanahan, Irish Lab June 2017.
## <danieldmcclanahan@gmail.com>

library(flowCore, warn.conflicts = F)
library(CytobankAPI, quietly = T, warn.conflicts = F)
library(flowWorkspace, quietly = T)
library(flowUtils)
library(CytoML, warn.conflicts = F)
library(tools)
library(digest)
library(Rtsne, warn.conflicts = F)
library(emdist, warn.conflicts = F)
library(spade, quietly = T, warn.conflicts = F, verbose = F)
library(gdata, warn.conflicts = F)
library(magrittr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(ggplot2)


### Read/write different representations of flow data.

## NOTE: we only grab a single dataset from each input file, even if more
## exist. multiple datasets in an fcs file is deprecated and "implementors are
## being discouraged to do so." it may exist, but shouldn't be supported.

read_fcs_raw <- function (fname) {
    ## TODO: are these the standard transformations for flowCore?
    read.FCS(fname, transformation = NULL, truncate_max_range = F)
}

read_fcs <- function (fname) {
    flowFrame_obj <- read_fcs_raw(fname)
    fname %>% read_fcs_raw %>% exprs %>% as.data.frame
}

read_text_file <- function (fname, ...) {
    tryCatch(
        read.table(fname, header = header, sep = sep, ...),
        error = function (e) {
            read.table(fname, skip = 1, ...)
        }
    )
}

read_txt <- function (fname, ...) {
    read_text_file(fname, header = T, sep = "\t", ...)
}

read_csv <- function (fname, ...) {
    read_text_file(fname, header = T, sep = ",", ...)
}

## If you need to do something weird when reading files, set this to your read
## function.
read_file_fun <- NULL
## Modify this list to recognize more filetypes. Filetypes are case-insensitive.
read_file_funs_map <- list(fcs = read_fcs,
                           csv = read_csv,
                           txt = read_txt)
read_file <- function (fname) {
    ## TODO: does any kind of data cleaning make sense here? see ../README.md
    ## TODO: consider having a cache for this function if files are reused a lot
    ext <- file_ext(fname)
    ext_lowercase <- tolower(ext)
    reader_fun <-
        if (!is.null(read_file_fun)) {
            match.fun(read_file_fun)
        } else if (ext_lowercase %in% names(read_file_funs_map)) {
            match.fun(read_file_funs_map[[ext_lowercase]])
        } else {
            stop(sprintf("unrecognized extension '%s' for file '%s'",
                         ext, fname))
        }
    reader_fun(fname)
}

write_txt <- function (frame, fname) {
    write.table(frame, sep = '\t', row.names = F, file = fname)
}

write_csv <- function (frame, fname) {
    write.table(frame, sep = ',', row.names = F, file = fname)
}

## TODO: check validity of data in df? against params/description?
make_flowFrame <- function (exprs, parameters, description) {
    exprs <- as.matrix(exprs)
    args <- list(exprs = quote(exprs))
    if (hasArg(parameters)) {
        args <- c(args, list(parameters = quote(parameters)))
    }
    if (hasArg(description)) {
        args <- c(args, list(description = quote(description)))
    }
    do.call("flowFrame", args, envir = environment())
}

write_flowFrame <- function (flow_obj, fname) {
    write.FCS(flow_obj, fname)
}

get_fcs_data_files <- function (dir = getwd(),
                                 pattern = "\\.fcs$",
                                 recursive = F) {
    list.files(path = dir, pattern = pattern, all.files = T, full.names = T,
               recursive = recursive, no.. = T)
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

## TODO: source this transform and why it's useful
asinh_transform <- function (x) {
    asinh(x / 5)
}

## TODO: make config file to canonicalize column names!
## Return data frame which contains only marker data.
## sample_n: equal sampling
## sample_frac: proportional sampling
fcs_data_cols <- function (
    frame,
    excl_pats = list(numeric_pat, sne_pat),
    excl_preds = list(is_all_integer)
    ) {
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



### Manipulate fcs data.

shared_markers <- function (frames) {
    ## TODO: find less common columns and check if they're mistakes
    ## TODO: if columns are close but not the same (e.g. levenshtein), show a
    ## warning
    frames %>% lapply(colnames) %>% Reduce(f = intersect)
}

do_tsne <- function (infiles, n,
                     markers = NULL, transform_with = asinh_transform,
                     verbose = T, ...) {
    num_f <- length(infiles)
    suffix <- sprintf("_visne_events_%s.fcs", n)
    outfiles <- gsub("\\.fcs$", suffix, infiles, perl = T)
    if (verbose) {
        cat(sprintf("sampling %s rows per file for %s files...\n",
                    n, num_f))
    }
    sampled_frames <- infiles %>% lapply(function (file) {
        if (verbose) {
            cat(sprintf("reading %s...\n", file))
        }
        read_file(file) %>% sample_n(size = n)
    })
    on_markers <-
        if (!is.null(markers)) {
            markers
        } else {
            if (verbose) {
                cat("guessing markers to join on...\n")
            }
            sampled_frames %>% lapply(fcs_data_cols) %>% shared_markers
        }
    if (verbose) {
        cat(sprintf("joining on markers:\n[%s]\n",
                    paste0(on_markers, collapse = ", ")))
        cat("performing visne...\n")
    }
    ## join into one big table to perform tsne
    joined <- sampled_frames %>% lapply(function (df) {
        ## select the shared columns and squash the raw data
        ## throws error if a frame doesn't have all the markers specified
        df[,on_markers] %>% mutate_all(transform_with)
    }) %>% Reduce(f = rbind)
    tsne_frame <- Rtsne(
        joined, dims = 2, check_duplicates = F, verbose = verbose, ...)
    if (verbose) {
        cat("writing visne results to files...\n")
    }
    ## write to fcs files and check that they exist
    lapply(1:num_f, function (i) {
        sampled_frame <- sampled_frames[[i]]
        outfile <- outfiles[i]
        if (verbose) {
            cat(sprintf("writing %s (%s/%s)...\n", outfile, i, num_f))
        }
        ## find the tsne axes that belong to us
        tsne_start <- ((i - 1) * n) + 1
        tsne_range <- tsne_start:(tsne_start + n - 1)
        ## add tsne axes and convert to flowFrame
        sampled_frame %>%
            mutate(tSNE1 = tsne_frame$Y[tsne_range,1],
                   tSNE2 = tsne_frame$Y[tsne_range,2]) %>%
            make_flowFrame %>% write.FCS(filename = outfile)
    })
    stopifnot(all(file.exists(outfiles)))
    outfiles
}



### Analyze hierarchies of populations in a dataset.

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

sort_by_component <- function (strs, split_by,
                               orders = list(), split_fixed = T,
                               value = T) {
    splits <- if (split_fixed) {
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

## TODO: parallelize this!
emd_fcs <- function (files, outfile,
                     max_iterations = 10,
                     verbose = T) {
    stopifnot(!any(duplicated(files)))
    n <- length(files)
    if (verbose) {
        cat(sprintf("reading in %s files...\n", n))
    }
    mats <- lapply(files, function (file) {
        read_file(file) %>% select(c(tSNE1, tSNE2)) %>% as.matrix
    })
    output <- matrix(vector(mode = "double", length = n * n), nrow = n)
    if (verbose) {
        cat(sprintf("starting pairwise emd on %s files...\n", n))
    }
    for (i in 1:n) {
        if (verbose) {
            cat(sprintf("row %s/%s\n", i, n))
        }
        i_mat <- mats[[i]]
        i_rows <- dim(i_mat)[1]
        if (i > 1) {
            for (j in 1:(i - 1)) {
                if (verbose) {
                    cat(sprintf("column %s/%s\n", j, n))
                }
                output[i,j] <- output[j,i]
            }
        }
        output[i,i] <- 0
        if (verbose) {
            cat(sprintf("column %s/%s\n", i, n))
        }
        if (i < n) {
            for (j in (i + 1):n) {
                if (verbose) {
                    cat(sprintf("column %s/%s\n", j, n))
                }
                j_mat <- mats[[j]]
                j_rows <- dim(j_mat)[1]
                output[i,j] <- emdw(i_mat, rep(1, i_rows),
                                    j_mat, rep(1, j_rows),
                                    max.iter = max_iterations)
            }
        }
    }
    colnames(output) <- files
    rownames(output) <- files
    write.table(output, outfile, sep = ",")
    stopifnot(file.exists(outfile))
    outfile
}

## look at xtabs/ftable/table() and summary/summarize/aggregate/group_by()



### Pull data from cytobank.

## TRUE unless file exists, is readable, and has the right size and contents
invalid_fcs_dl <- function (fcs_info) {
    fcs_info %$%
        {
            (file.exists(filename) &
             file.access(filename, mode = 4) == 0 &
             fileSize == file.size(filename) &
             md5sum == md5sum(filename))
        } %>%
        ## any NA or F means file is inaccessible or corrupted
        ## (e.g. by an interrupted download)
        { is.na(.) | !(.) }
}

## return fcs filenames
download_all_fcs <- function (session, exp_id, verbose = T) {
    fcs_info <- fcs_files.list(session, exp_id) %>%
        ## fcs_files.list returns a data frame with list columns; undo that
        mutate_all(unlist) %>%
        select(id, filename, md5sum, fileSize) %>%
        mutate(humansize = humanReadable(fileSize, width = 4))
    to_dl <- fcs_info[invalid_fcs_dl(fcs_info),]
    size_all <- sum(to_dl$fileSize)
    size_cur <- 0
    if (verbose) {
        cat(sprintf("total (unzipped) size of files to download: %s\n",
                    humanReadable(size_all, width = 4)))
    }
    num_to_dl <- dim(to_dl)[1]
    if (num_to_dl != 0) {
        for (i in 1:num_to_dl) {
            row <- to_dl[i,]
            if (verbose) {
                cat(sprintf("%s downloaded (%s%% done)\n",
                            humanReadable(size_cur, width = 4),
                            round(size_cur / size_all * 100, digits = 1)))
                cat(sprintf("downloading %s (unzipped size: %s, id: %s)...\n",
                            row$filename, row$humansize, row$id))
            }
            taken <- system.time(
                fcs_files.download_zip(session, exp_id, row$id) %>% unzip)
            cat(sprintf("download of %s took %s seconds\n",
                        row$filename, round(taken[3], 2)))
            if (any(invalid_fcs_dl(row))) {
                stop(sprintf(paste("download of %s is corrupt",
                                   "-- check your internet connection"),
                             row$filename))
            }
            size_cur <- size_cur + row$fileSize
        }
    }
    fcs_info$filename
}

download_fcs_cytobank <- function (experiments, verbose = T) {
    if (length(experiments) == 0) { return(character()) }
    library(getPass)
    session <- authenticate(readline("site (<site>.cytobank.org): "),
                            readline("username: "),
                            getPass(forcemask = T))
    Reduce(f = c, init = character(),
           x = lapply(experiments, function (exp_id) {
               download_all_fcs(session, exp_id, verbose = verbose)
           }))
}

download_gates <- function (session, exp_id) {
    gates.gatingML_download(session, exp_id)
}

apply_gates <- function (gates_xml, fcs_files) {
    cytobank2GatingSet(gates_xml, fcs_files)
}



### Comparative clustering.

daniel_read_fcs <- function (fname) {
    read.FCS(fname, transformation = NULL, truncate_max_range = F) %>%
        exprs %>% as.data.frame %>%
        mutate(manual_gates = as.integer(manual_gates))
}

num_term_pops <- function (frame) {
    pops <- frame$manual_gates %>% unique %>% sort
    k <- length(pops)
    stopifnot(all(pops == 1:k))
    k
}

f_measure_clusters <- function (frame) {
    cls <- frame$cluster %>% unique %>% sort
    k <- length(cls)
    stopifnot(all(cls == 1:k))
    lapply(cls, function (cl_ind) {
        mode_pop <- frame[frame$cluster == cl_ind,] %>% .$manual_gates %>%
            table %>% sort(decreasing = T) %>%
            names %>% as.integer %>% .[1]
        n_relevant <- sum(frame$manual_gates == mode_pop &
                          frame$cluster == cl_ind)
        precision <- n_relevant / sum(frame$cluster == cl_ind)
        recall <- n_relevant / sum(frame$manual_gates == mode_pop)
        2 * precision * recall / (precision + recall)
    }) %>% unlist %>% mean
}

cluster_kmeans <- function (frame, k) {
    clustering <- frame %>% select(c(tSNE1, tSNE2)) %>% kmeans(centers = k)
    frame %>% mutate(cluster = clustering$clust)
}

cluster_spade <- function (files) {
    k <- lapply(files, read_file) %>% Reduce(f = bind_rows) %>% num_term_pops
    cat(sprintf("k = %s\n", k))
    suppressWarnings(
        SPADE.driver(files = files, cluster_cols = c("tSNE1", "tSNE2"),
                     transforms = NULL, k = k,
                     downsampling_target_percent = .01))
    clustered <- sprintf("%s.density.fcs.cluster.fcs", files)
    lapply(clustered, read_fcs) %>% Reduce(f = bind_rows)
}
