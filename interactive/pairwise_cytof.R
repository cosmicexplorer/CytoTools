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
emd_tsne_outfile <- "emd_tsne.csv"
emd_heatmap_outfile <- "heatmap_emd.pdf"
mem_outfile <- "mem.csv"
mem_rmsd_outfile <- "mem_rmsd.csv"
mem_heatmap_outfile <- "heatmap_mem_rmsd.pdf"

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
write.csv(emd_pairwise_matrix, emd_tsne_outfile)

pdf(emd_heatmap_outfile)
## read the EMD comparison matrix back from emd_outfile. check.names = FALSE
## ensures it doesn't mess with our column names's spaces or anything.
emd_matrix_fromfile <- as.matrix(read.csv(emd_tsne_outfile, row.names = 1,
                                          check.names = FALSE))
CytoTools::plot_pairwise_comparison(emd_matrix_fromfile, col = color_palette)
dev.off()



### MEM analysis
pheno_data <- CytoTools::normalize_pheno_channels_dataset(
    cytof_data, ref = CytoTools::read_cyto_file("./iPSCs.fcs"),
    transform_fun = CytoTools::asinh_transform,
    channel_name_ops = c(
        CytoTools::non_pheno_channel_name_patterns,
        list("^NA$" = NA,
             "^IR$" = NA,
             "^CISPLATIN$" = NA)))

markers <- pheno_data$shared_channels
message(sprintf("joining on %s markers: [%s]",
                length(markers), paste(markers, collapse = ", ")))
mem_df <- CytoTools::calc_mem(
    pheno_data$pop_list, pheno_data$ref,
    IQRthresh = 0.5, scale_limit = 10)

## write this to file so we save all our hard work
write.csv(mem_df, mem_outfile)
## create distance matrix (with euclidean distance) from scaled MEM vectors
write.csv(as.matrix(dist(mem_df)), mem_rmsd_outfile)

pdf(mem_heatmap_outfile)
## read the MEM RMSD comparison matrix back from file so we check that we wrote
## it correctly -- exact same as using as.matrix(dist(mem_df))
mem_rmsd_fromfile <- as.matrix(
    read.csv(mem_rmsd_outfile, row.names = 1, check.names = FALSE))
CytoTools::plot_pairwise_comparison(
    mem_rmsd_fromfile, with_dendrograms = TRUE,
    ## play with the number of colors in this variable (at the top of this file)
    ## if the coloration seems weird or the color key is blank
    col = color_palette,
    ## play with these to adjust axis label sizes
    cexRow = .25, cexCol = .25, margin = c(5, 5))
dev.off()
