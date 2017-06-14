library(magrittr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(flowCore, warn.conflicts = F)
library(CytobankAPI, quietly = T, warn.conflicts = F)
library(flowWorkspace, quietly = T)
library(flowUtils)
library(CytoML)
library(tools)


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

## site="irishlab", username=<>, password=<>
cytobank_request <- function (..., req) {
    cy_sesh <- authenticate(...)
    tryCatch({
        match.fun(req)(cy_sesh)
    },
    finally = authentication.logout(cy_sesh))
}

fcs_successful_downloads <- function (fcs_tbl) {
    with(fcs_tbl, (
        file.exists(path) &
        file.access(path, mode = 4) &
        fileSize == file.size(path) &
        md5sum == md5sum(path)))
}

fetch_fcs_unzip <- function (..., fcs_tbl, dir = getwd()) {
    failed <- fcs_tbl %>% fcs_failed_downloads
    fcs_tbl %>%
        mutate(path = {

        }
                   (fileSize == file.size(path) &
                    md5sum == md5sum(path) &
                    ))
        {
            fcs_failed_downloads %<>%
                mutate(path = {
                    if (n() > 0) {
                        fcs_zip_path <- fcs_files.download_zip(
                            ..., directory = dir, fcs_files = id) %>%
                            as.vector(mode = "character")
                        unzip(fcs_zip_path, files = filename, exdir = dir)
                    }
                }) %T>% (
                    fcs_failed_downloads %>% summarize(stopifnot(n() == 0)))
        }
}

download_all_fcs <- function (auth, exp_id, dir = getwd()) {
    fcs_files_info <- fcs_files.list(auth, exp_id) %>%
        select(id, filename, md5sum, fileSize) %>%
        mutate(path = file.path(dir, filename)) %>%
        ## ensure it's a data frame with vectors
        as.data.frame %>% mutate_all(unlist) %>%

        group_by(fetched =
                     fileSize == file.size(path) &
                     md5sum == md5sum(path)) %>%
        mutate()

    if (length(fcs_to_fetch) != 0) {
        fcs_files_zip <- fcs_files.download_zip(
            auth, experiment_id = exp_id, fcs_files = id, directory = dir) %>%
            as.vector(mode = "character")
        stopifnot(length(fcs_files_zip) == 1)

        unzip(fcs_files_zip, )
    }







    ## throw unless we can download everything
    ## return downloaded file paths
}

dl_check_fcs <- function (auth, exp_id, fcs_id, file_md5) {
    fcs_dl <- fcs_files.download(
        auth, experiment_id = exp_id, fcs_file_id = fcs_id)

}

get_gates_pops_set <- function (auth, exp_id, dir = getwd(),
                                without_fcs_dl = FALSE) {
    ## download every fcs file from the experiment
    if (!without_fcs_dl) {
        all_fcs_zipped <- fcs_files.download_zip(auth, experiment_id = exp_id)
        unzip(all)
    }
    ## download gatingml as xml
    gating_file <- gates.gatingML_download(auth, experiment_id = exp_id)
    cytobank2GatingSet(gating_file, list.files(path = "./5-2"))
}
