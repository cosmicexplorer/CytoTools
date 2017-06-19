library(flowCore, warn.conflicts = F)
library(CytobankAPI, quietly = T, warn.conflicts = F)
library(flowWorkspace, quietly = T)
library(flowUtils)
library(CytoML)
library(gdata, warn.conflicts = F)
library(magrittr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)


### Read/write different representations of flow data.

## NOTE: we only grab a single dataset from each input file, even if more
## exist. multiple datasets in an fcs file is deprecated and "implementors are
## being discouraged to do so." it may exist, but shouldn't be supported.

read_fcs <- function(fname) {
    flowFrame_obj <- read.FCS(
        ## TODO: are these the standard transformations for flowCore?
        fname, transformation = NULL, truncate_max_range = F)
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

write_flowFrame <- function (frame, fname) {
    write.FCS(frame, fname)
}



### Clean fcs data.

numeric_pat <- "^[[:digit:]+]$"
sne_pat <- "[sS][nN][eE]"

is_all_integer <- function (vec) {
    ## vec is a column of a data frame produced by read_clean_fcs
    ## read.FCS only produces double columns
    stopifnot(is.vector(vec) && is.double(vec))
    ## this benchmarked as the fastest method
    all(vec == as.integer(vec))
}

## TODO: make name_changes (able to) read from config file!
clean_frame <- function (name_changes = list(),
                         excl_pats = list(numeric_pat, sne_pat),
                         excl_preds = list(is_all_integer),
                         frame) {
    canonized <- frame %>% colnames %>%
        Reduce(x = name_changes, init = .,
               f = function (cols, fun) { lapply(cols, fun) }) %>%
        set_colnames(x = frame, value = .)

    names_filtered <- canonized %>% colnames %>%
        Reduce(x = excl_pats, init = .,
               f = function (names, pat) {
                   matched <- grepl(pat, names, perl = T)
                   names[matched] <- NA
                   names
               }) %>%
        set_colnames(x = canonized, value = .)


    col_value_filtered <- names_filtered %>% select_if({
        for (pred in excl_preds) {
            ## short-circuit failure
            if (pred(.)) { return(FALSE) }
        }
        TRUE
    })

    col_value_filtered
}



### Manipulate fcs data.

shared_markers <- function (frames) {
    ## TODO: find less common columns and check if they're mistakes
    ## TODO: if columns are close but not the same (e.g. levenshtein), show a
    ## warning
    frames %>% colnames %>% Reduce(f = intersect, x = .)
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
