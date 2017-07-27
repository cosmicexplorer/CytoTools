### Compute pairwise Earth Mover's Distance (EMD) and MEM RSMD and generate
### heatmaps.
## Written by Danny McClanahan <danieldmcclanahan@gmail.com> with the Irish Lab
## at the Vanderbilt University Medical Center Department of Cancer Biology.

library(CytoTools)

### Configure
## `color_palette`: Color palette for heatmap.
color_palette <- colorRampPalette(
    c("#24658C", "#5CBAA7", "#9ED2A4", "#E2E998", "#FBF7BF", "#FDDC86",
      "#F8A05A", "#EF6342", "#D43E4F")
)(n = 50)
## `emd_outfile`: filename for output csv with EMD pairwise comparisons
emd_outfile <- "emd_out.csv"
## `emd_heatmap_outfile`: output file containing EMD heatmap as a pdf
emd_heatmap_outfile <- "heatmap_emd.pdf"
## `mem_outfile`: filename for output csv with MEM RMSD pairwise comparisons
mem_outfile <- "mem_out.csv"
## `mem_heatmap_outfile`: output file containing MEM heatmap as a pdf
mem_heatmap_outfile <- "heatmap_mem.pdf"

data_files <- CytoTools::fcs_file_paths(path = ".", pattern = "\\.fcs$")

data_files_sorted <- CytoTools::sort_by_component(
    data_files,
    split_by = "_",
    orders = list(c(), c("pre", "3wk", "12wk", "6m")))

## Read files into a list of data frames.
cytof_data <- lapply(data_files_sorted, CytoTools::read_cyto_file)
## Names are used to mark each file's row and column in the final
## heatmap. The below uses ~regex magic~ to make the names look a little
## prettier for a specific format of filename -- just do:
##
## > names(cytof_data) <- data_files_sorted
##
## for the easiest option.
names(cytof_data) <- CytoTools::replace_matches(
    data_files_sorted,
    list("^[\\./]*" = "",
         "\\.fcs$" = "",
         "viSNE" = "",
         "^_+|_+$" = "",
         "_+" = " "))



### EMD analysis on viSNE axes
tsne_matrices <- lapply(cytof_data, function (fcs_df) {
    as.matrix(fcs_df[,c("tSNE1", "tSNE2")])
})

emd_pairwise_matrix <- CytoTools::pairwise_emd(tsne_matrices)

## write the EMD comparison matrix to emd_outfile
write.csv(emd_pairwise_matrix, emd_outfile)

pdf(emd_heatmap_outfile)
## read the EMD comparison matrix back from emd_outfile. check.names = FALSE
## ensures it doesn't mess with our column names's spaces or anything.
emd_matrix_fromfile <- as.matrix(read.csv(emd_outfile, row.names = 1,
                                          check.names = FALSE))
CytoTools::plot_pairwise_comparison(emd_matrix_fromfile, col = color_palette)
dev.off()



### MEM analysis
mem_pairwise_matrix <- CytoTools::pairwise_mem_rmsd(
    cytof_data, ref_pop = CytoTools::read_cyto_file("./iPSCs.fcs"),
    ## by default, channel names are transformed with ?toupper when calling
    ## ?CytoTools::normalize_channels (which is where the ... arguments to
    ## ?CytoTools::pairwise_mem_rmsd go to), so patterns are written in
    ## uppercase.
    channel_name_ops = c(
        CytoTools::non_pheno_channel_name_patterns,
        list("^NA$" = NA,
             "^IR$" = NA,
             "^CISPLATIN$" = NA)))

## write the MEM RMSD comparison matrix to mem_outfile
write.csv(mem_pairwise_matrix, mem_outfile)

pdf(mem_heatmap_outfile)
## read the MEM RMSD comparison matrix back from mem_outfile
mem_matrix_fromfile <- as.matrix(read.csv(mem_outfile, row.names = 1,
                                          check.names = FALSE))
CytoTools::plot_pairwise_comparison(
    mem_matrix_fromfile,
    col = color_palette,
    cexRow = .25, cexCol = .25, margin = c(5, 5),
    dendro = TRUE)
dev.off()
