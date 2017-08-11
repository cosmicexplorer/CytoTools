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
emd_tsne_means_outfile <- "emd_tsne_means.csv"
emd_tsne_variances_outfile <- "emd_tsne_variances.csv"
emd_tsne_medians_outfile <- "emd_tsne_medians.csv"
emd_heatmap_outfile <- "heatmap_emd.pdf"
mem_outfile <- "mem.csv"
mem_rmsd_outfile <- "mem_rmsd.csv"
mem_heatmap_outfile <- "heatmap_mem_rmsd.pdf"

sorted_cytof_data <- CytoTools::read_files_sorted(
    path = ".",
    extensions = c("fcs"),
    split_by = c("_", " "),
    orders = list(c(), c("pre", "3wk", "12wk", "6m")))


### EMD analysis on viSNE axes
tsne_matrices <- lapply(sorted_cytof_data, function (fcs_df) {
    as.matrix(fcs_df[,c("tSNE1", "tSNE2")])
})

emd_pairwise_analysis <- CytoTools::pairwise_emd(
    tsne_matrices, downsample_rows = 200L, comparison_runs = 10L,
    summary_funs = list(
        mean = mean,
        variance = var,
        median = median),
    verbose_timing = TRUE)

## save the EMD summary stats matrices to file so we don't lose our work
write.csv(emd_pairwise_analysis$mean, emd_tsne_means_outfile)
write.csv(emd_pairwise_analysis$variance, emd_tsne_variances_outfile)
write.csv(emd_pairwise_analysis$median, emd_tsne_medians_outfile)

pdf(emd_heatmap_outfile)
CytoTools::plot_pairwise_comparison(
    emd_pairwise_analysis$mean,
    pcnt_similarity_scale = TRUE,
    ## play with the number of colors in this variable (at the top of this file)
    ## if the coloration seems weird or the color key is blank
    col = color_palette)
dev.off()

## you can quickly see a histogram to see if any variances seem too high
hist(emd_pairwise_analysis$variance / emd_pairwise_analysis$mean)

## or, you can plot a heatmap and see which pairs of files have the greatest
## variance. this isn't super useful -- if the variance is too high on any of
## them, you should increase the number of rows sampled per file or increase the
## number of comparison runs
CytoTools::plot_pairwise_comparison(
    emd_pairwise_analysis$variance,
    col = color_palette)


### MEM analysis

pheno_channel_names <- c(
    "CCR4", "CCR5", "CCR7", "CD16", "CD19", "CD20", "CD27",
    "CD28", "CD38",
    "CD44", "CD45", "CD45RA", "CD45RO", "CD56",
    "CD69", "CD8", "CXCR3", "CXCR5", "HLA-DR", "ICOS", "PD-1", "TCRgd",
    "TIM3")
pheno_data <- CytoTools::normalize_pheno_channels_dataset(
    pheno_channel_names, sorted_cytof_data,
    ref = CytoTools::read_cyto_file("./iPSCs.fcs"),
    transform_fun = asinh_transform(5))

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
    mem_rmsd_fromfile, with_dendrograms = TRUE, col = color_palette,
    ## play with these to adjust axis label sizes
    cexRow = .25, cexCol = .25, margin = c(5, 5))
dev.off()
