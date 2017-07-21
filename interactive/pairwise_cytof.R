### Compute pairwise Earth Mover's Distance (EMD) and MEM RSMD and generate
### heatmaps.
## Written by Danny McClanahan, Irish Lab June 2017.
## <danieldmcclanahan@gmail.com>

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
## Set working directory, if you need to.
setwd(".")

data_files <- list.files(pattern = "\\.fcs$", ignore.case = TRUE,
                         all.files = TRUE, full.names = TRUE, recursive = FALSE,
                         no.. = TRUE)

data_files_sorted <- CytoTools::sort_by_component(
    data_files,
    split_by = "_",
    orders = list(c(), c("pre", "3wk", "12wk", "6m")))

name_repls <- c("^[cC][dD]([0-9]+)\\-[0-9]+$" = "CD\\1")
cyto_data_frames <- CytoTools::process_cyto_dataset(
    data_files_sorted, name_repls)

CytoTools::pairwise_emd(cyto_data_frames, emd_outfile, max_iterations = 10)

pdf(emd_heatmap_outfile)
CytoTools::plot_pairwise_comparison(emd_outfile, color_palette = color_palette)
dev.off()

mem_ref_pop <- NULL
CytoTools::pairwise_mem_rmsd(
    cyto_data_frames, mem_outfile, ref_pop = mem_ref_pop)

pdf(mem_heatmap_outfile)
CytoTools::plot_pairwise_comparison(
    mem_outfile, color_palette = color_palette, dendro = T)
dev.off()
