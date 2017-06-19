library(flowCore, warn.conflicts = F)
library(CytobankAPI, quietly = T, warn.conflicts = F)
library(flowWorkspace, quietly = T)
library(flowUtils)
library(CytoML)
library(tools)
library(Rtsne)
library(spade, quietly = T, warn.conflicts = F)
library(gdata, warn.conflicts = F)
library(magrittr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)


### Read/write different representations of flow data.

## NOTE: we only grab a single dataset from each input file, even if more
## exist. multiple datasets in an fcs file is deprecated and "implementors are
## being discouraged to do so." it may exist, but shouldn't be supported.

read_fcs_raw <- function (fname) {
    ## TODO: are these the standard transformations for flowCore?
    read.FCS(fname, transformation = NULL, truncate_max_range = F)
}

read_fcs <- function(fname) {
    flowFrame_obj <- read_fcs_raw(fname)
    ## TODO: does any kind of data cleaning make sense here? see ../README.md
    as.data.frame(exprs(flowFrame_obj))
}

read_csv <- function (fname, header = TRUE, sep = ',', ...) {
    read.table(fname, header = header, sep = sep, ...)
}

read_txt <- function (fname, header = TRUE, sep = '\t', ...) {
    read.table(fname, header = header, sep = sep, ...)
}

read_file <- function (fname) {
    ext <- file_ext(fname)
    switch(ext,
           fcs = read_fcs(fname),
           csv = read_csv(fname),
           txt = read_txt(fname),
           stop(sprintf("unrecognized file extension '%s' for file '%s'",
                        ext, fname)))
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

do_tsne <- function (files, n,
                     transform = asinh_transform,
                     outdir = getwd(), verbose = T,
                     ...) {
    outfiles <- gsub(
        "\\.fcs$", sprintf("_sampled_%s.fcs", n), basename(files), perl = T)
    sampled <- files %>% setNames(., .) %>% lapply(function (file) {
        read_file(file) %>% sample_n(size = n)
    })
    ## get data columns shared between all
    on_markers <- sampled %>% lapply(fcs_data_cols) %>% shared_markers
    tsne <- sampled %>% lapply(function (df) {
        ## select the shared columns and squash the raw data
        df[,on_markers] %>% mutate_all(match.fun(transform))
    }) %>% Reduce(f = rbind) %>%
        ## join into one big table to perform tsne
        Rtsne(check_duplicates = F, verbose = verbose, ...)
    print(tsne)
    with_sne <- lapply(1:length(sampled), function (i) {
        ## find the tsne axes that belong to us
        tsne_start <- ((i - 1) * n) + 1
        tsne_range <- tsne_start:(tsne_start + n - 1)
        ## add tsne axes and convert to flowFrame
        sampled[[i]] %>%
            mutate(tSNE1 = tsne$Y[tsne_range,1],
                   tSNE2 = tsne$Y[tsne_range,2]) %>%
            make_flowFrame
    })
    ## write to fcs files and check that they exist
    with_sne %>% as("flowSet") %>%
        ## write.flowSet is still "experimental"
        {
            suppressWarnings(
                write.flowSet(x = ., outdir = outdir, filename = outfiles))
        }
    stopifnot(all(file.exists(outfiles)))
    outfiles
}



### Analyze hierarchies of populations in a dataset.

## look at xtabs/ftable/table() and summary/summarize/aggregate/group_by()


### Pull data from cytobank.

## check that file exists, is readable, and has the right size and contents
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
        cat(sprintf("total size of files to download: %s\n",
                    humanReadable(size_all, width = 4)))
    }
    for (i in 1:dim(to_dl)[1]) {
        row <- to_dl[i,]
        if (verbose) {
            cat(sprintf("%s downloaded (%s%% done)\n",
                        humanReadable(size_cur, width = 4),
                        round(size_cur / size_all * 100, digits = 1)))
            cat(sprintf("downloading %s (zipped size: %s, id: %s)...\n",
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
    fcs_info$filename
}

download_gates <- function (session, exp_id) {
    gates.gatingML_download(session, exp_id)
}

apply_gates_fcs <- function (gates_xml, fcs_files) {
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

cluster_spade <- function (files, k) {
    SPADE.driver(files = files, cluster_cols = c("tSNE1", "tSNE2"),
                 transforms = NULL, k = k)
    clustered <- sprintf("%s.density.fcs.cluster.fcs", files)
    lapply(clustered, read_fcs) %>% Reduce(f = bind_rows)
}
