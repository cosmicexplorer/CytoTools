### Run viSNE on a set of files.
## Written by Danny McClanahan, Irish Lab June 2017.
## <danieldmcclanahan@gmail.com>

source("./fcs.R")

### Get data files to run analysis on.
## `cytobank_experiments`: List of experiment IDs on cytobank.
##   Download all fcs files from all experiments in this list to the current
##   directory, and store the file names in `cytobank_fcs`. You don't need this
##   if you already have your data files.
##
##   If `length(cytobank_experiments) > 0`, the terminal will prompt you for
##   your cytobank URL, username, and password. The password prompt loads the
##   library "getPass", which will ensure it is not echoed to the terminal. You
##   may need to call `install.packages("getPass")`.
cytobank_experiments <- list()
cytobank_fcs <- download_fcs_cytobank(cytobank_experiments, verbose = verbose)

## data_files: Char vector of files to read data from.
##   Files can be binary fcs files or text files with headers. ".txt" files will
##   be read as TSV (sep = "\t"), and ".csv" files as CSV (sep = ",").
##
##   If the file is a text file (.txt or .csv) and starts with a blank line,
##   this script will recognize that and skip the initial blank line.
##
##   You can manually set data_files as well. For example:
##   > data_files <- c("file1.fcs", "file2.fcs")
data_files <- list.files(pattern = "\\.fcs$", ignore.case = TRUE,
                         all.files = TRUE, full.names = TRUE, recursive = FALSE,
                         no.. = TRUE)

## with_tsne: Char vector of files produced by do_tsne().
##   do_tsne() takes a random sample of n events per input file (without
##   replacement), performs a viSNE analysis on the markers specified, and
##   returns the generated filenames with tSNE axes tSNE1 and tSNE2.
##
##   (Note that viSNE is performed locally for now, so it may require some wait
##   time to run if n is large.)
##
##   If markers = NULL, do_tsne() will guess what the markers to use
##   should be. It will display what those markers are.
##
##   transform_with specifies how to transform the data before performing
##   viSNE. The default is to use the asinh(x / 5) transformation. Set
##   transform_with = identity if your data is already transformed.
##
##   Other arguments are passed with ... to Rtsne. See ?Rtsne for other
##   parameters of the viSNE analysis.
with_tsne <- do_tsne(data_files, n = 10000, markers = NULL,
                     perplexity = 30, theta = 0.5, max_iter = 1000)

## with_tsne now contains the files with tSNE axes added!
