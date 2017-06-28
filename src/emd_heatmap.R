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
## `emd_outfile`: filename for output csv with EMD pairwise comparisons
emd_outfile <- "emd_out.csv"
## `heatmap_outfile`: output file containing EMD heatmap as a pdf
heatmap_outfile <- "heatmap_fcs_out.pdf"
## Set working directory, if you need to.
setwd(".")

## data_files: Char vector of files to read data from.
##   Files can be binary fcs files or text files with headers. ".txt" files will
##   be read as TSV (sep = "\t"), and ".csv" files as CSV (sep = ",").
##
##   NOTE: Files should have two tSNE axes labeled tSNE1 and tSNE2!!!
##
##   If the file is a text file (.txt or .csv) and starts with a blank line,
##   this script will recognize that and skip the initial blank line.
##
##   You can manually set data_files as well. For example:
##   > data_files <- c("file1.fcs", "file2.fcs")
data_files <- list.files(pattern = "\\.fcs$", ignore.case = TRUE,
                         all.files = TRUE, full.names = TRUE, recursive = FALSE,
                         no.. = TRUE)

## sort_files_by_component(): splits filenames into pieces and sorts them.
##   sort_files_by_component() uses the string in `split_by` to break filenames
##   into components.
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
##   > sort_files_by_component(
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
data_files_sorted <- sort_files_by_component(
    data_files,
    split_by = "",
    orders = list())

## emd_fcs(): Run pairwise EMD on input files and produce CSV.
##   The resuts are stored in `emd_outfile`.
##
##   `max_iterations` is the number of iterations to perform when computing EMD.
##   Increasing this value *typically* does not change the result at all.
emd_fcs(data_files_sorted, emd_outfile, max_iterations = 10)

## emd_outfile has row names in column 1
emd_matrix <- as.matrix(read.csv(row.names = 1))

pdf(heatmap_outfile)
heatmap.2(emd_matrix, Rowv = F, Colv = F, dendrogram = "none",
          col = color_palette, trace = "none")
dev.off()
