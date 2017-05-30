library(flowCore)
library(magrittr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)

with_connection <- function (conn, expr, name = "con", pf = parent.frame()) {
    sub_conn <- substitute(conn)
    sub_expr <- substitute(expr)
    toDo <- bquote(with(setNames(list(.(sub_conn)), .(name)),
                        .(sub_expr)))
    eval(toDo, pf)
}

## see readFCSheader in IO.R in ~/tools/flowCore/R
## use with_connection(file("marrow_Ungated_viSNE.fcs", "r"), readFCSheader(con))

read_clean_fcs <- function(fname) {
    flowFrame_obj <- read.FCS(
        ## TODO: are these the standard transformations for flowCore?
        fname, transformation = NULL, truncate_max_range = F)
    ## TODO: does any kind of data cleaning make sense here? see ../README.md
    as.data.frame(exprs(flowFrame_obj))
}

## read list of fcs filenames into data frames
read_fcs_files <- function(fnames) {
    names_list <- as.list(fnames)
    lapply(names_list, read_clean_fcs)
}

## TODO: convert binary fcs to text and back!

squish_expression <- Vectorize(function(x) asinh(x / 5))
inv_squish_expression <- Vectorize(function(y) 5 * sinh(y))

add_squished <- function (frame, channels,
                          sep = "S-", squish = squish_expression) {
    squished_names <- paste(sep, colnames(frame)[which(channels)], sep = "")
    squished <- squish(frame[,channels])
    colnames(squished) <- squished_names
    cbind(frame[,!channels], frame[,channels], squished)
}

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

filter_columns <- function (name_changes, excl_pats, excl_preds, df) {
    canonized <- df %>% colnames %>%
        Reduce(x = name_changes, init = .,
               f = function (cols, fun) { lapply(cols, fun) }) %>%
        set_colnames(x = df, value = .)

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

process_intersect_fcs <- function (frames,
                                   ## TODO: add config file for this!
                                   name_changes = list(),
                                   excl_pats = list(numeric_pat, sne_pat),
                                   excl_preds = list(is_all_integer)) {
    filtered <- lapply(frames, function (df) {
        filter_columns(name_changes, excl_pats, excl_preds, df)
    })

    ## TODO: find less common columns and check if they're mistakes
    ## TODO: if columns are close but not the same (levenshtein), show a warning
    common_markers <- filtered %>% colnames %>% Reduce(f = intersect, x = .)
    lapply(filtered, function (df) df[,common_markers])
}
