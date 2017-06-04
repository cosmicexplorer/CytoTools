library(flowCore)
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

## TODO: convert binary fcs to text and back!
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
    ## this was found to be the fastest way to do this with benchmarking
    ## (not like that matters here)
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
