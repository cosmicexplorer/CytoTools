source("fcs.R")
library(stringr)
library(stringdist)
library(magrittr, warn.conflicts = F)

check_invalid_changes <- function(orig_cols, new_cols) {
    if (length(new_cols) != length(orig_cols)) {
        stop(sprintf(paste(
            "more or less canonicalized columns than input channels:\n",
            "fcs columns (%d): [%s]\n", "canonicalized(%d): [%s]", sep = ""),
            length(orig_cols), str(orig_cols),
            length(new_cols), str(new_cols)))
    } else if (0 != anyDuplicated(new_cols)) {
        stop(sprintf(paste(
            "some channels were converted to the same canonical column",
            "[input->canon]:\n%s\n"),
            paste(orig_cols, new_cols, sep = "->", collapse = "\n")))
    }
}

## pop is a list with frame, name, and ref_level, as in pop_ref_match().
## policy() returns the canonicalized column name, or NA.
canonicalize_channels <- function(pop, policy) {
    fr <- pop$frame                     # original frame
    orig_cols <- colnames(fr)

    ## we canonicalize the column names of the data frame and record all changes
    ## made in changes, where the key is the canonicalized name and the value is
    ## the previous name
    ## do a right reduce so no rev required (this is NOT recursive)
    canon_to_prev <- Reduce(
        x = orig_cols,
        init = list(cols = list(),
                    changes = list(),
                    noncanonized = list()),
        right = T,
        f = function(channel, acc) {
            canon <- policy(channel)
            if (!is.na(canon)) {
                acc$changes[[canon]] = channel
                new <- canon
            } else {
                acc <- prepend_to_key_if(acc, 'noncanonized', channel)
                new <- channel
            }
            list(cols = c(new, acc$cols),
                 changes = acc$changes,
                 noncanonized = acc$noncanonized)
        })

    ## halt if two columns canonicalized to same value
    new_cols <- canon_to_prev$cols
    check_invalid_changes(orig_cols, new_cols)

    pop$frame <- setNames(fr, new_cols)
    pop$changes <- canon_to_prev$changes
    pop$noncanonized <- canon_to_prev$noncanonized

    pop
}

## does NOT mutate lst
prepend_to_key_if <- function(lst, key, val) {
    if (is.list(lst[[key]])) {
        lst[[key]] <- c(val, lst[[key]])
        lst
    } else { NULL }
}

quantile_levels <- c(.25, .5, .75)

## pops as in pop_ref_match(), but with "changes" field from
## canonicalize_channels()
add_pop_stats <- function(pops) {
    all_events <- sum(
        vapply(pops, function(p) dim(p$frame)[1], vector("integer", 1)))
    ## adds "abundance" and "quantiles" fields
    lapply(pops, function(pop) {
        pop$abundance <- dim(pop$frame)[1] / all_events
        pop$quantiles <- apply(pop$frame, 2, function (channel) {
            ## TODO: is this the correct definition of quartiles/median?
            ## "quantile" has a "type" argument which produces different output
            quantile(channel, quantile_levels)
        })
        pop
    })
}

## *_quartiles should be vector of length 3 with .25/.5/.75 quantiles in order,
## or the first quartile, the median, then the third quartile.
## MEM is defined as in "Characterizing cell subsets using marker enrichment
## modeling", Nature Methods, Diggins et al. 2017
## TODO: use the MEM function from the paper's code instead of redoing it here?
MEM_from_quartiles <- function (pop_quartiles, ref_quartiles) {
    stopifnot(length(pop_quartiles) == 3 && length(ref_quartiles) == 3)

    IQR_POP <- pop_quartiles[3] - pop_quartiles[1]
    MAG_POP <- pop_quartiles[2]
    IQR_REF <- ref_quartiles[3] - ref_quartiles[1]
    MAG_REF <- ref_quartiles[2]

    sign <- if (MAG_POP - MAG_REF < 0) { -1 } else { 1 }
    MEM_base <- abs(MAG_POP - MAG_REF) + (IQR_REF / IQR_POP) - 1
    MEM_base * sign
}

MEM_row_names <- c("Q1", "Median", "Q3", "MEM")

## compare a single pop to all the reference populations.
## assume pop names are unique -> remove elem with pop$name from refs if exists.
## adds "comparisons" field to pop, which contains a list of all ref pops, but
## with "quantiles" field restricted ONLY to common channels, as well as:
## - MEM -> [data.frame] MEM scores for the common columns
compare_pop_to_refs <- function(pop, refs) {
    refs_without_cur <- Filter(x = refs, function (ref) ref$name != pop$name)
    pop$comparisons <- lapply(refs_without_cur, function (ref) {
        common_cols <- intersect(colnames(pop$frame), colnames(ref$frame))
        MEM_ref_to_pop <- sapply(common_cols, function (col) {
            pop_quartiles <- pop$quantiles[,col]
            ref_quartiles <- ref$quantiles[,col]
            MEM_score <- MEM_from_quartiles(pop_quartiles, ref_quartiles)
            c(ref_quartiles, MEM_score)
        })
        ref$MEM <- as.data.frame(row.names = MEM_row_names, MEM_ref_to_pop)
        ref
    })
    pop
}

## policy is as in canonicalize_channels()
## pops is a list of lists, where each sublist has the following keys:
## - frame -> [data.frame] table containing results of cytometry for population
## - filename -> [string] file path to fcs file containing population
## - name -> [string] name of population
## - ref_level -> [string] type of population (see below)
pop_ref_match <- function (processed) {
    pops <- processed$pops

    ## adds "abundance" and "quantiles" fields
    pops_with_stats <- add_pop_stats(pops)

    ## split into three classes:
    ## - no_ref -> population is not compared to others as a ref pop
    ## - ref_global -> population compared as ref, not added to MEM heatmap
    ## - ref_relative -> population compared as ref, added to MEM heatmap (with
    ##   MEM label of all zeroes)
    init <- list(no_ref = list(),
                 ref_global = list(),
                 ref_relative = list())
    ## do a right reduce so no rev required (this is NOT recursive)
    pop_classes <- Reduce(
        x = pops_with_stats,
        init = init,
        right = T,
        f = function (pop, acc) {
            ## use list(pop) because otherwise list is flattened
            res <- prepend_to_key_if(acc, pop$ref_level, list(pop))
            stopifnot(!is.null(res))
            res
        })

    pops_non_ref <- pop_classes$no_ref
    pops_global_ref <- pop_classes$ref_global
    pops_rel_ref <- pop_classes$ref_relative

    ## TODO: generate heatmaps!
    all_pops <- c(pops_non_ref, pops_global_ref, pops_rel_ref)
    all_refs <- c(pops_global_ref, pops_rel_ref)
    processed$analyzed <- lapply(all_pops, function (pop) {
        compare_pop_to_refs(pop, all_refs)
    })

    processed
}

select_matches_from <- function(input, reg_list, use_perl, case_sensitive) {
    cols <- if (case_sensitive) { input } else { toupper(input) }
    pats <- if (case_sensitive) { reg_list } else { toupper(reg_list) }
    Reduce(
        x = pats,
        init = rep(F, times = length(cols)),
        f = function(matches, reg) {
            matches | grepl(reg, cols, fixed = !use_perl, perl = use_perl)
        })
}

## removes columns of only integers, as well as any columns with names
## corresponding to any element of remove_column_regexps (a list of strings)
strip_integer_and_matched_cols <- function (frame, filename,
                                            remove_column_regexps,
                                            strip_ints, use_perl, case,
                                            verbose) {
    orig_cols <- colnames(frame)

    ## get columns which don't match regexps
    matched_by_name <- select_matches_from(
        orig_cols, remove_column_regexps, use_perl, case)
    name_cols_removed <- orig_cols[matched_by_name]
    write_verbose(verbose,
                  sprintf("columns removed by regexp (in file '%s'): [%s]",
                          filename, paste(name_cols_removed, collapse = ", ")))
    by_name <- which(!matched_by_name)

    if (!strip_ints) {
        list(data = frame[,by_name],
             name_cols_removed = name_cols_removed,
             int_cols_removed = c())
    } else {
        ## get columns which aren't all integers
        stripped <- Reduce(
            x = by_name, init = list(keep = list(), discard = list()),
            right = T, f = function (col, acc) {
                if (!all(frame[,col] %% 1 == 0)) {
                    prepend_to_key_if(acc, 'keep', col)
                } else {
                    prepend_to_key_if(acc, 'discard', col)
                }
            })
        int_cols_removed <- orig_cols[unlist(stripped$discard)]
        write_verbose(verbose,
                      sprintf("columns of all integers removed: [%s]",
                              paste(int_cols_removed, collapse = ", ")))
        list(data = frame[,unlist(stripped$keep)],
             name_cols_removed = name_cols_removed,
             int_cols_removed = int_cols_removed)
    }
}

default_canonicalize_policy <- function(x) { NA }

compare_with_case <- function (s1, s2, case_sensitive) {
    if (case_sensitive) {
        s1 == s2
    } else {
        toupper(s1) == toupper(s2)
    }
}

apply_when <- function (val, bool, fn) {
    if (bool) { fn(val) } else { val }
}

set_cond <- function (bool, yes_val, no_val) {
    if (bool) { yes_val } else { no_val }
}

file_canonicalize_policy <- function(infiles, case, in_canon) {
    infiles_fmt <- paste(infiles, collapse = ", ")
    channel_matches <- list()
    canon <- apply_when(in_canon, !case, function (c) lapply(c, toupper))
    case_str <- set_cond(case, "on", "off")
    for (infile in infiles) {
        matches <- str_match(read_lines_device(infile), "([^:]+):([^:]+)")
        for (rowNum in 1:dim(matches)[1]) {
            row <- matches[rowNum,]
            raw_channel_name <- row[2]
            channelName <- apply_when(raw_channel_name, !case, toupper)
            raw_canonical <- row[3]
            canonical <- apply_when(raw_canonical, !case, toupper)
            already_matched_canonical <- channel_matches[[channelName]]
            if (!is.null(already_matched_canonical)) {
                if (already_matched_canonical == canonical) {
                    warning(sprintf(
                        paste("the channel '%s' was mapped (in file '%s')",
                              "to the canonical name '%s' more than once",
                              "(with case sensitivity turned %s)",
                              "in files [%s]"),
                        raw_channel_name, infile, canonical, case_str,
                        infiles_fmt))
                } else {
                    stop(sprintf(
                        paste("the channel '%s' was mapped (in file '%s')",
                              "to the canonical name '%s',",
                              "which is different than the previous '%s'",
                              "(with case sensitivity turned %s)",
                              "in files [%s]"),
                        raw_channel_name, infile, canonical,
                        already_matched_canonical, case_str, infiles_fmt))
                }
            } else {
                if (is.na(match(canonical, canon))) {
                    same_canon <- canon[ain(canon, canonical, maxDist = 2)]
                    if (length(same_canon) != 0) {
                        warning(sprintf(paste(
                            "the canonical name '%s' for channel '%s'",
                            "(in file '%s') is similar to these other",
                            "canonical names",
                            "(with case sensitivity turned %s):",
                            "[%s]\nthis may indicate an",
                            "error in files [%s]"),
                            raw_canonical, raw_channel_name, infile, case_str,
                            paste(same_canon, collapse = "', '"), infiles_fmt))
                    }
                    canon <- c(canonical, canon)
                }
                channel_matches[[channelName]] = canonical
            }
        }
    }

    list(canon = canon,
         policy = function(raw_channel_name) {
             channelName <- apply_when(raw_channel_name, !case, toupper)
             mapped <- channel_matches[[channelName]]
             if (is.null(mapped)) { NA } else { mapped }
         })
}

read_lines_device <- function (file) {
    handle <- file(file, raw = T)
    open(handle, open = "r")
    res <- readLines(handle)
    close(handle)
    res
}

write_lines_device <- function (inp, file) {
    handle <- file(file, raw = T)
    open(handle, open = "w")
    writeLines(inp, handle)
    close(handle)
    NULL
}

join_column_regexp_files <- function (files) {
    if (is.null(files)) { return(c()) }
    lapply(files, read_lines_device) %>% unlist
}

merge_canons <- function (canons, case) {
    Reduce(
        x = apply_when(canons, !case, toupper),
        init = c(),
        f = function (acc, cur) {
            if (cur %in% acc) { acc } else { c(cur, acc) }
        })
}

## TODO: add "cluster index" as nullable string key for elements of pops_fcs;
## split rows of the fcs file into pops based on the cluster index column

## pops_fcs is a list of lists, where each sublist has the following keys:
## - filename -> [string] path to fcs file
## - pop_name -> [string (nullable)] name of population in output
## - ref_level -> [string] see pop_ref_match()
## FIXME: read multiple data segments in a single fcs file to separate pops!
## channel_canonicalize_file -> [string] path to file with canonize mappings
## remove_column_regexps are used with the grep() function, which accepts
## partial matches -- be aware of this!
process_fcs <- function (pops_fcs, in_canon, canonize_files,
                         remove_column_regexps_files, strip_ints, use_perl,
                         case, verbose, clean, change_names) {
    canonization <- file_canonicalize_policy(canonize_files, case, in_canon)
    policy <- canonization$policy
    remove_column_regexps <- join_column_regexp_files(
        remove_column_regexps_files)
    pops <- lapply(pops_fcs, function (pop) {
        ## use pop_name if given, otherwise use filename as pop_name
        pop_name <- if (!is.null(pop$pop_name)) { pop$pop_name }
                    else { basename(pop$filename) }
        stripped <- strip_integer_and_matched_cols(
            read_clean_fcs(pop$filename, clean, change_names),
            pop$filename, remove_column_regexps, strip_ints, use_perl, case,
            verbose)
        list(frame = stripped$data,
             filename = pop$filename,
             name = pop_name,
             ref_level = pop$ref_level,
             name_cols_removed = stripped$name_cols_removed,
             int_cols_removed = stripped$int_cols_removed)
    })

    ## ensure names are unique
    if (0 != anyDuplicated(lapply(pops, function (p) p$name))) {
        stop(sprintf("names of input pops must be unique:\n%s",
                     str(pops)))
    }

    ## modifies colnames of frame and adds "changes" key to each sublist
    canonicalized <- lapply(pops, function (p) canonicalize_channels(p, policy))

    full_canon <- merge_canons(
        unlist(list(canonization$canon,
                    lapply(canonicalized, function (p) p$noncanonized))),
        case)

    list(pops = canonicalized,
         canon = full_canon,
         remove_column_regexps_files = remove_column_regexps_files,
         canonize_files = canonize_files)
}

format_changes <- function(changes, filename = NULL) {
    lapply(names(changes), function (n) {
        fname_str <- if (is.null(filename)) { "" }
                     else { sprintf(" (in '%s')", filename) }
        sprintf("%s -> %s%s", changes[[n]], n, fname_str)
    })
}

## apply to output of "pop_ref_match" to print output
## can use sink() to divert output to file
print_text_output <- function(analysis) {
    pop_ref_matches <- analysis$analyzed
    for (pop in pop_ref_matches) {
        cat(sprintf("%s (abundance=%.2f%%, type=%s, file=%s):\n",
                    pop$name, pop$abundance * 100, pop$ref_level,
                    pop$filename))
        print(pop$quantiles)
        for (ref in pop$comparisons) {
            cat(sprintf("\tref: %s (abundance=%.2f%%, type=%s, file=%s):\n",
                        ref$name, ref$abundance * 100, ref$ref_level,
                        ref$filename))
            print(ref$MEM)
        }
        cat("-------------------\n")
    }
    all_removed_regexp <- pop_ref_matches %>% lapply(function (pop) {
        lapply(pop$name_cols_removed, function (col) {
            sprintf("%s (in file '%s')", col, pop$filename)
        })}) %>% unlist
    cat(sprintf("all removed by regexp (from files [%s]):\n%s\n",
                paste(analysis$remove_column_regexps_files, collapse = ", "),
                paste(all_removed_regexp, collapse = "\n")))

    all_removed_int <- pop_ref_matches %>% lapply(function (pop) {
        lapply(pop$int_cols_removed, function (col) {
            sprintf("%s (in file '%s')", col, pop$filename)
        })}) %>% unlist
    cat(sprintf("all integer columns removed:\n%s\n",
                paste(all_removed_int, collapse = "\n")))

    all_changes_list <- Reduce(
        x = pop_ref_matches, init = list(), function (acc, pop) {
            c(acc, format_changes(pop$changes, pop$filename))
        })
    cat(sprintf("all changes to canonical (from files [%s]):\n%s\n",
                paste(analysis$canonize_files, collapse = ", "),
                paste(all_changes_list, collapse = "\n")))
}

write_verbose <- function (verb, str) {
    if (verb) {
        write(str, stderr())
    }
}

perform_verbose <- function (msg, verb, expr) {
    ev <- substitute(expr)
    if (!verb) { eval.parent(ev) }
    else {
        time <- system.time(res <- eval.parent(ev))
        write(sprintf("(%s took %f seconds)", msg, time[3]), stderr())
        res
    }
}
