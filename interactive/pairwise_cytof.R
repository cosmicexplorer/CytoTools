### Compute pairwise Earth Mover's Distance (EMD) and MEM RSMD and generate
### heatmaps.
## Written by Danny McClanahan, Irish Lab June 2017.
## <danieldmcclanahan@gmail.com>

library(CytoTools)
library(magrittr)

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
## Set working directory, if you need to.
setwd(".")

data_files <- CytoTools::fcs_file_paths(path = ".")

data_files_sorted <- CytoTools::sort_by_component(
    data_files,
    split_by = "_",
    orders = list(c(), c("pre", "3wk", "12wk", "6m")))

name_repls <- c("^[cC][dD]([0-9]+)\\-[0-9]+$" = "CD\\1")
cytof_data <- lapply(data_files_sorted, function (file_path) {
    CytoTools::read_cyto_file(file_path, rx_replace = name_repls)
}) %>% setNames(data_files_sorted)

emd_pairwise_matrix <- CytoTools::pairwise_emd(cytof_data)
## write the EMD comparison matrix to emd_outfile
write.table(emd_pairwise_matrix, emd_outfile, sep = ',')

pdf(emd_heatmap_outfile)
## read the EMD comparison matrix back from emd_outfile
emd_matrix_fromfile <- emd_outfile %>% read.csv(row.names = 1) %>% as.matrix
CytoTools::plot_pairwise_comparison(
    emd_matrix_fromfile, color_palette = color_palette)
dev.off()

mem_ref_pop <- NULL
marker_names <- CytoTools::get_marker_names(
    c(mem_ref_pop, cytof_data),
    name_filter = c(CytoTools::non_pheno_channel_name_patterns,
                    ## matches intercalator, e.g. "NA-191"
                    "^NA-[0-9]+$",
                    ## matches cisplatin, e.g. "Cisplatin-195"
                    "^Cisplatin-[0-9]+$"))
mem_pairwise_matrix <- CytoTools::pairwise_mem_rmsd(
    cytof_data, markers = marker_names, ref_pop = mem_ref_pop)
## write the MEM RMSD comparison matrix to mem_outfile
write.table(mem_pairwise_matrix, mem_outfile, sep = ',')

pdf(mem_heatmap_outfile)
## read the MEM RMSD comparison matrix back from mem_outfile
mem_matrix_fromfile <- mem_outfile %>% read.csv(row.names = 1) %>% as.matrix
CytoTools::plot_pairwise_comparison(
    mem_matrix_fromfile, color_palette = color_palette, dendro = TRUE)
dev.off()
