library(flowCore, warn.conflicts = F)
library(dplyr, warn.conflicts = F)

read_fcs <- function (fname) {
    read.FCS(fname, transformation = NULL, truncate_max_range = F) %>%
        exprs %>% as.data.frame %>%
        mutate(manual_gates = as.integer(manual_gates))
}

num_term_pops <- function (frame) {
    pops <- frame$manual_gates %>% unique %>% sort
    k <- length(pops)
    stopifnot(all(pops == 1:k))
    k
}

f_measure_clusters <- function (frame, clusters) {
    stopifnot(length(clusters) == dim(frame)[1])
    cls <- clusters %>% unique %>% sort
    k <- length(cls)
    stopifnot(all(cls == 1:k))
    lapply(cls, function (cl_ind) {
        clust <- frame[clusters == cl_ind,]
        mode_pop <- table(clust$manual_gates) %>% names %>% as.integer %>% .[1]
        relevant <- clust[clust$manual_gates == mode_pop,]
        precision <- relevant / dim(clust)[1]
        in_mode_pop <- frame[frame$manual_gates == mode_pop,]
        recall <- relevant / dim(in_mode_pop)[1]
        2 * precision * recall / (precision + recall)
    }) %>% unlist %>% mean
}

cluster_kmeans <- function (frame, k) {
    frame %>% select(c(tSNE1, tSNE2)) %>% kmeans(centers = k)
}
