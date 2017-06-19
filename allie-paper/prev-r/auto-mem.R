#!/usr/bin/env Rscript

suppressMessages(library(argparse))

source("pop-ref-match.R")

arg_parser <- ArgumentParser(
    description = 'Process FCS files for MEM scores',
    epilog = paste(
        '-c and -r can be specified multiple times to use multiple files.',
        'The script reads each line of column match files as a pattern. By',
        'default, these are interpretened as fixed strings to be',
        'partially matched, ignoring case. Adding -I stops ignoring case, and',
        '-P makes the script interpret these as PCRE regular expressions.',
        'See http://www.regular-expressions.info/tutorial.html for a tutorial,',
        'or https://regex101.com/ to test regular expressions.',
        collapse = ' '))
arg_parser$add_argument(
               '-v', '--verbose', action = 'store_true',
               help = 'print timing and debug output')
arg_parser$add_argument(
               '-f', '--column-matches-file', action = 'append',
               help = 'remove columns with names matching patterns in file')
arg_parser$add_argument(
               '-s', '--no-strip-ints', action = 'store_false',
               help = 'do NOT strip integer-only columns from input files')
arg_parser$add_argument(
               '-P', '--use-perl-regexps', action = 'store_true',
               help = 'use PCRE regular expressions for patterns')
arg_parser$add_argument(
               '-I', '--case-sensitive', action = 'store_true',
               help = 'make patterns case sensitive')
arg_parser$add_argument(
               '-r', '--ref', action = 'append',
               help = 'use file as a reference population')
arg_parser$add_argument(
               '-c', '--canonize-mapping-file', action = 'append',
               help = 'read channel canonicalization from this file')
arg_parser$add_argument(
               '-C', '--channel-canon-input-file', default = NULL,
               help = 'file with one canonical channel name per line')
arg_parser$add_argument(
               '-O', '--channel-canon-output-file', default = NULL,
               help = paste('where to store all canonical channels from this',
                            'analysis, in addition to the -C argument',
                            '(if supplied)'))
arg_parser$add_argument(
               'files', metavar = 'FILE...',
               help = 'non-ref FCS files to be analyzed')

analyze <- function(arguments) {
    infiles <- lapply(arguments$files, function (f) {
        list(filename = f, ref_level = "no_ref")
    })
    ref_files <- lapply(arguments$ref, function (r) {
        list(filename = r, ref_level = "ref_global")
    })
    input <- append(infiles, ref_files)

    canon_inf <- arguments$channel_canon_input_file
    in_canon <- if (is.null(canon_inf)) { list() }
                else { as.list(read_lines_device(canon_inf)) }

    write_verbose(arguments$verbose, "reading/processing fcs files...")
    pops <- perform_verbose(
        "processing", arguments$verbose,
        process_fcs(
            input, in_canon, arguments$canonize_mapping_file,
            arguments$column_matches_file, arguments$no_strip_ints,
            arguments$use_perl_regexps, arguments$case_sensitive,
            arguments$verbose, T, F))

    canon_outf <- arguments$channel_canon_output_file
    if (!is.null(canon_outf)) {
        write_lines_device(unlist(pops$canon), canon_outf)
    }

    write_verbose(arguments$verbose, "performing analysis...")
    perform_verbose("analysis", arguments$verbose, pop_ref_match(pops))
}

if (!interactive()) {
    args <- arg_parser$parse_args()
    output <- analyze(args)
    print_text_output(output)
}

## can load this and call analyze function directly
## example command line:
## ./auto-mem.R -v -f example-regexps-file -f example-reg-file-2 \
##              -c example-canon-file -c example-canon-file-2 \
##              -r '1305007 D0 PB.fcs' -r '1305007 D8.fcs' \
##              '1305007 D0 BM.fcs'
## corresponds to the call:
## print_text_output(analyze(
##     column_matches_file = list("example-regexps-file", "example-reg-file-2"),
##     canonize_mapping = list("example-canon-file", "example-reg-file-2"),
##     files = list("1305007 D0 BM.fcs"),
##     ref = list("1305007 D0 PB.fcs", "1305007 D8.fcs"),
##     verbose = T, no_strip_ints = T, use_perl_regexps = F,
##     case_sensitive = F, channel_canon_input_file = NULL,
##     channel_canon_output_file = NULL))
## can call library(profvis) and wrap the above in profvis() for profiling
