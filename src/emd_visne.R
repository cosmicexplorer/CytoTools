### Run viSNE on a set of files, compute pairwise Earth Mover's Distance (EMD),
### and generate a heatmap.
source("./fcs.R")

### Configure
## `color_palette`: Color palette for heatmap.
color_palette <- colorRampPalette(
    c("#24658C", "#5CBAA7", "#9ED2A4", "#E2E998", "#FBF7BF", "#FDDC86",
      "#F8A05A", "#EF6342", "#D43E4F")
)(n = 299)
## `heatmap_outfile`: output file containing EMD heatmap as a pdf
heatmap_outfile <- "heatmap_fcs_out.pdf"
## `verbose`: Whether to show progress for each call.
verbose <- T
## Set working directory, if you need to.
setwd(".")

### Get data files to run analysis on.
## `cytobank_experiments`: List of experiment IDs on cytobank.
##   Download all fcs files from all experiments in this list to the current
##   directory, and store the file names in `cytobank_fcs`. Don't touch this if
##   you already have your data files.
##
##   If this is non-empty, the terminal will prompt you for your cytobank URL,
##   username, and password. The password prompt loads the library "getPass",
##   which will ensure it is not echoed to the terminal.
cytobank_experiments <- list()
cytobank_fcs <- download_fcs_cytobank(cytobank_experiments, verbose = verbose)

## `data_files`: Char vector of files to read data from.
##   Files can be binary fcs files or text files with headers. ".txt" files will
##   be read as TSV (sep = "\t"), and ".csv" files as CSV (sep = ",").
##
##   If the file is a text file and starts with a blank line, this script will
##   recognize that and skip the initial blank line.
data_files <- list.files(pattern = "\\.fcs$", ignore.case = T,
                         all.files = T, full.names = T, recursive = F,
                         no.. = T)[1:13]

## `with_tsne`: Char vector of files produced by do_tsne().
##   do_tsne() takes a random sample of `n` events per input file (without
##   replacement) and performs a viSNE analysis on the markers specified. See
##   ?Rtsne for other parameters of the viSNE analysis. The resulting output
##   files will have axes tSNE1 and tSNE2.
##
##   (Note that viSNE is performed locally for now, so it may require some wait
##   time to run if `n` is large.)
##
##   If `markers` is NULL, do_tsne() will guess what the markers to use
##   should be. If `verbose` is TRUE, it will display what those markers are.
##
##   `transform` specifies how to transform the data before
##   performing viSNE. `asinh_transform` performs the asinh(x / 5)
##   transformation. Use `identity` if your data is already transformed.
##
##   do_tsne() will reuse its output if the input files' content and `n` have
##   not changed. It checks file content changing with md5sum(). Set
##   use_existing = F to turn this off (e.g. if you change a viSNE parameter).
with_tsne <- do_tsne(data_files, markers = NULL,
                     n = 1000, perplexity = 10, theta = 0.5, max_iter = 200,
                     transform = asinh_transform,
                     use_existing = T, verbose = verbose)

## `pairwise_emd_table`: File containing pairwise EMD matrix.
##   `max_iterations` is the number of iterations to perform when computing EMD.
##   Increasing this value *typically* does not change the result at all.
##
##   `use_existing` is the same as in do_tsne().
pairwise_emd_table <- emd_fcs(with_tsne, max_iterations = 10,
                              use_existing = T, verbose = verbose)

emd_matrix <- as.matrix(read.csv(pairwise_emd_table, row.names = 1))
## TODO: sort here!

pdf(heatmap_outfile)
heatmap(emd_matrix, Rowv = NA, Colv = NA, col = color_palette)
dev.off()
