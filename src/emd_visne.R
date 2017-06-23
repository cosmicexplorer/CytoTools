### Run viSNE on a set of files, compute pairwise Earth Mover's Distance (EMD),
### and generate a heatmap.
## Written by Danny McClanahan, Irish Lab June 2017.
## <danieldmcclanahan@gmail.com>

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
verbose <- TRUE
## Set working directory, if you need to.
setwd(".")

### Get data files to run analysis on.
## `cytobank_experiments`: List of experiment IDs on cytobank.
##   Download all fcs files from all experiments in this list to the current
##   directory, and store the file names in `cytobank_fcs`. Don't touch this if
##   you already have your data files.
##
##   If `length(cytobank_experiments) > 0`, the terminal will prompt you for
##   your cytobank URL, username, and password. The password prompt loads the
##   library "getPass", which will ensure it is not echoed to the terminal. You
##   may need to call `install.packages("getPass")`.
cytobank_experiments <- list()
cytobank_fcs <- download_fcs_cytobank(cytobank_experiments, verbose = verbose)

## `data_files`: Char vector of files to read data from.
##   Files can be binary fcs files or text files with headers. ".txt" files will
##   be read as TSV (sep = "\t"), and ".csv" files as CSV (sep = ",").
##
##   If the file is a text file (.txt or .csv) and starts with a blank line,
##   this script will recognize that and skip the initial blank line.
data_files <- list.files(pattern = "\\.fcs$", ignore.case = TRUE,
                         all.files = TRUE, full.names = TRUE, recursive = FALSE,
                         no.. = TRUE)

## `sort_by_component()`: splits filenames into pieces and sorts them.
##   sort_by_component() uses the string in `split_by` to break filenames into
##   components.
##
##   Splitting "A:B:C" by ":" returns c("A", "B", "C")). We say that the string
##   "A:B:C" has components "A" at index 1, "B" at index 2, and "C" at 3.
##
##   `orders` is a list of character vectors. If a component shows up here at
##   the correct index, its file is ordered before other files. Among files
##   which have matching components at each index, the order is determined by
##   the order in the character vector at that index.
##
##   Example:
##
##   > data_files <- c("MB004_6m_panel2.fcs", "MB004_3wk_panel2.fcs",
##                     "MB004_12wk_panel2.fcs", "MB004_pre_panel2.fcs",
##                     "MB005_12wk_panel2.fcs")
##   > sort_by_component(
##       data_files,
##       split_by = "_",
##       orders = list(c("MB005"), c("pre", "3wk", "12wk", "6m")))
##   [1] "MB005_12wk_panel2.fcs" "MB004_pre_panel2.fcs"  "MB004_3wk_panel2.fcs"
##   [4] "MB004_12wk_panel2.fcs" "MB004_6m_panel2.fcs"
##
##   "MB005" sorts before "MB004" in the first component above because "MB005"
##   is mentioned in `orders` at index 1, but "MB004" isn't. If we did it again
##   with `orders` = list(c(), c("pre", "3wk", "12wk", "6m")), the first
##   component would be sorted alphabetically, so "MB005" would come after
##   "MB004".
##
##   With `split_by` = "" and `orders` = list(), `data_files` is simply sorted
##   alphabetically by file name.
data_files <- sort_by_component(
    data_files,
    split_by = "",
    orders = list())

data_files <- sort_by_component(
    data_files,
    split_by = "_",
    orders = list(c(), c("pre", "3wk", "5wk", "6wk", "12wk", "6m")))

## `with_tsne`: Char vector of files produced by do_tsne().
##   do_tsne() takes a random sample of `n` events per input file (without
##   replacement) and performs a viSNE analysis on the markers specified. See
##   ?Rtsne for other parameters of the viSNE analysis. The resulting output
##   files will have axes tSNE1 and tSNE2.
##
##   (Note that viSNE is performed locally for now, so it may require some wait
##   time to run if `n` is large.)
##
##   If `markers` = NULL, do_tsne() will guess what the markers to use
##   should be. If `verbose` is TRUE, it will display what those markers are.
##
##   `transform` specifies how to transform the data before
##   performing viSNE. `asinh_transform` performs the asinh(x / 5)
##   transformation. Use `identity` if your data is already transformed.
##
##   do_tsne() will reuse its output if the input files' content and `n` have
##   not changed. It checks file content changing with md5sum(). Set
##   `use_existing` = FALSE to turn this off. For example if you change a viSNE
##   parameter, such as perplexity, do_tsne() will use the old output if the
##   input files and `n` are the same, which is NOT what you want. So you can
##   set `use_existing` = FALSE to have it perform viSNE with the new
##   perplexity. If `verbose` = TRUE, do_tsne() will say whether it is reusing
##   an old output.
with_tsne <- do_tsne(data_files, n = 10000, markers = NULL,
                     perplexity = 30, theta = 0.5, max_iter = 1000,
                     transform = asinh_transform,
                     use_existing = TRUE, verbose = verbose)

## `pairwise_emd_table`: File containing pairwise EMD matrix.
##   `max_iterations` is the number of iterations to perform when computing EMD.
##   Increasing this value *typically* does not change the result at all.
##
##   `use_existing` is the same as in do_tsne(). If `use_existing` is TRUE and
##   the input files are the same, but in a different order, emd_fcs() will not
##   recompute the matrix, but it will rearrange the rows and columns according
##   to the new order.
pairwise_emd_table <- emd_fcs(with_tsne,
                              max_iterations = 10,
                              use_existing = TRUE, verbose = verbose)

emd_matrix <- as.matrix(read.csv(pairwise_emd_table, row.names = 1))

pdf(heatmap_outfile)
heatmap(emd_matrix, Rowv = NA, Colv = NA, col = color_palette)
dev.off()
