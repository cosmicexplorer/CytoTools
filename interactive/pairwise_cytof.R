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

emd_pairwise_analysis <- CytoTools::pairwise_emd(
    tsne_matrices, downsample_rows = 100L, comparison_runs = 5L,
    summary_funs = list(mean = mean, variance = var),
    verbose_timing = TRUE)

x <- CytoTools::pairwise_emd(
    tsne_matrices, downsample_rows = 50L, comparison_runs = 2L,
    summary_funs = list(mean = mean, variance = var))

## save the EMD mean and variance matrices to file so we don't lose our work
write.csv(emd_pairwise_analysis$mean, emd_tsne_means_outfile)
write.csv(emd_pairwise_analysis$variance, emd_tsne_variances_outfile)

pdf(emd_heatmap_outfile)
## read the EMD comparison matrix back from file. check.names = FALSE
## ensures it doesn't remove spaces or other characters in our column names

emd_means_fromfile <- as.matrix(read.csv(emd_tsne_means_outfile, row.names = 1,
                                         check.names = FALSE))

CytoTools::plot_pairwise_comparison(emd_means_fromfile, col = color_palette)

dev.off()

emd_variance_fromfile <- as.matrix(read.csv(emd_tsne_variances_outfile,
                                            row.names = 1, check.names = FALSE))

## you can quickly see a histogram to see if any variances seem too high
hist(emd_variance_fromfile)

## or, you can plot a heatmap and see which pairs of files have the greatest
## variance. this isn't super useful -- if the variance is too high on any of
## them, you should increase the number of rows sampled per file or increase the
## number of comparison runs
CytoTools::plot_pairwise_comparison(
    emd_variance_fromfile,
    ## don't do any scaling as with our distance matrices
    is_dist = FALSE,
    col = color_palette)


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
