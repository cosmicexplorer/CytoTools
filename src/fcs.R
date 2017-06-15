library(flowCore, warn.conflicts = F)
library(CytobankAPI, quietly = T, warn.conflicts = F)
library(flowWorkspace, quietly = T)
library(flowUtils)
library(CytoML)
library(curl)
library(foreach)
library(iterators)
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

## site="irishlab", username=<>, password=<>
cytobank_request <- function (..., req) {
    cy_sesh <- authenticate(...)
    tryCatch({
        match.fun(req)(cy_sesh)
    },
    finally = authentication.logout(cy_sesh))
}

invalid_fcs_dl <- function (fcs_list) {
    res <- with(fcs_list, (file.exists(path) &
                           file.access(path, mode = 4) &
                           fileSize == file.size(path) &
                           md5sum == md5sum(path)))
    is.na(res) | !res
}

get_fcs_zip <- function (..., fcs_info) {
    ids_str <- paste(fcs_info$id, collapse = ",")
    out <- sprintf("experiment_%s_fcs_%s.zip", exp_id, ids_str)
    url <- sprintf("%s/experiments/%s/fcs_files/download_zip?fcs_file_ids=%s",
                   session@site, exp_id, ids_str)
    h <- new_handle()
    handle_setheaders(
        h, "authorization" = sprintf("Bearer %s", session@auth_token))
    curl_download(url, out, quiet = F, handle = h)
    unzip(out) %T>%
        function (fcs_out) {
            stopifnot(setequal(normalizePath(fcs_out, mustWork = T),
                               fcs_info$path))
            stopifnot(!any(invalid_fcs_dl(fcs_info)))
        }
}

download_all_fcs <- function (..., directory = getwd()) {
    fcs_info <- fcs_files.list(...) %>%
        select(id, filename, md5sum, fileSize) %>%
        ## add where the files in experiment are expected to be
        mutate(path = {
            real_dir <- normalizePath(directory, mustWork = T)
            file.path(real_dir, filename)
        }) %>%
        ## ensure it's a data frame with vectors
        as.data.frame %>% mutate_all(unlist)

    foreach(i = iter(fcs_info, by = "row", chunksize = 2), .combine = "c") %:%
        when(invalid_fcs_dl(i)) %do%
        get_fcs_zip(..., fcs_files = i$id)
    message(sprintf("downloading the following:\n%s",
                    paste(fcs_info$path, collapse = "\n")))



    ## download POTENTIALLY HUGE zip file -- may fail
    fcs_zip_path <- fcs_files.download_zip(..., directory = directory) %>%
        as.vector(mode = "character")
    message(sprintf("zip file downloaded into '%s'", fcs_zip_path))

    ## get unzipped files
    unzip(fcs_zip_path, exdir = directory) %T>%
        ## throw unless we can downloaded / unzipped everything correctly
        ## return downloaded file paths
        function (fcs_out) {
            message(sprintf("unzipped files:\n%s",
                            paste(fcs_out, collapse = "\n")))
            stopifnot(setequal(normalizePath(fcs_out, mustWork = T),
                               fcs_info$path))

            successes <- with(fcs_info, (file.exists(path) &
                                         file.access(path, mode = 4) &
                                         fileSize == file.size(path) &
                                         md5sum == md5sum(path)))
            stopifnot(all(successes))
        }
}

get_gates_pops_set <- function (..., fcs_set = NULL) {
    ## download every fcs file from the experiment if not given
    if (is.null(fcs_set)) {
        fcs_set <- download_all_fcs(...)
    }
    ## download gatingml as xml
    gates.gatingML_download(...) %>%
        ## apply to specified fcs files
        cytobank2GatingSet(fcs_set)
}
