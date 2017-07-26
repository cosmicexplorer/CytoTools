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

cytof_data <- setNames(
    lapply(data_files_sorted, CytoTools::read_cyto_file),
    data_files_sorted)

tsne_matrices <- lapply(cytof_data, function (fcs_df) {
    as.matrix(fcs_df[,c("tSNE1", "tSNE2")])
})
emd_pairwise_matrix <- CytoTools::pairwise_emd(tsne_matrices)
## write the EMD comparison matrix to emd_outfile
write.table(emd_pairwise_matrix, emd_outfile, sep = ',')

pdf(emd_heatmap_outfile)
## read the EMD comparison matrix back from emd_outfile
emd_matrix_fromfile <- as.matrix(read.csv(emd_outfile, row.names = 1))
CytoTools::plot_pairwise_comparison(emd_matrix_fromfile, color_palette)
dev.off()

mem_ref_pop <- NULL
cytof_data_normalized <- CytoTools::normalize_channels(
    c(mem_ref_pop, cytof_data),
    channel_name_ops = c(CytoTools::non_pheno_channel_name_patterns,
                         "^NA$" = NA_character_,
                         "^CISPLATIN$" = NA_character_))
mem_pairwise_matrix <- CytoTools::pairwise_mem_rmsd(
    cytof_data_normalized, mem_ref_pop)
## write the MEM RMSD comparison matrix to mem_outfile
write.table(mem_pairwise_matrix, mem_outfile, sep = ',')

pdf(mem_heatmap_outfile)
## read the MEM RMSD comparison matrix back from mem_outfile
mem_matrix_fromfile <- as.matrix(read.csv(mem_outfile, row.names = 1))
CytoTools::plot_pairwise_comparison(mem_matrix_fromfile, color_palette,
                                    dendro = TRUE)
dev.off()
