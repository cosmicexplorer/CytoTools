### Routines for manipulating CyToF data and interacting with Cytobank.
## Written by Danny McClanahan, Irish Lab June 2017.
## <danieldmcclanahan@gmail.com>

library(flowCore, warn.conflicts = F)
library(CytobankAPI, quietly = T, warn.conflicts = F)
library(tools)
library(digest)
library(boot)
library(XML)
library(gplots)
library(hashmap)
library(Rtsne, warn.conflicts = F)
library(emdist, warn.conflicts = F)
library(spade, quietly = T, warn.conflicts = F, verbose = F)
library(gdata, warn.conflicts = F)
library(magrittr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)
library(ggplot2)
## local

### Read/write different representations of flow data.

## NOTE: we only grab a single dataset from each input file, even if more
## exist. Multiple datasets in an fcs file is deprecated and "implementors are
## being discouraged to do so." it may exist, but shouldn't be supported.

read_fcs_flowFrame <- function (fname) {
    read.FCS(fname, transformation = NULL, truncate_max_range = F)
}

read_fcs_cyto_frame <- function (fname) {
    fname %>% read_fcs_flowFrame %>% exprs %>% as.data.frame
}

read_text_cyto_frame <- function (fname, ...) {
    tryCatch(
        read.table(fname, header = T, ...),
        error = function (e) {
            read.table(fname, header = T, skip = 1, ...)
        }
    )
}

read_file <- function (fname) {
    ## TODO: does any kind of data cleaning make sense here? see ../README.md
    ## TODO: consider having a cache for this function if files are reused a lot
    ext <- file_ext(fname) %>% tolower
    switch(
        ext,
        fcs = read_fcs_cyto_frame(fname),
        csv = read_text_cyto_frame(fname, sep = ","),
        txt = read_text_cyto_frame(fname, sep = "\t"),
        stop(sprintf("unrecognized extension '%s' for file '%s'",
                     ext, fname)))
}

## TODO: check validity of data in df? against params/description?
make_flowFrame <- function (exprs, parameters, description) {
    exprs <- as.matrix(exprs)
    args <- list(exprs = quote(exprs))
    if (hasArg(parameters)) {
        args <- c(args, list(parameters = quote(parameters)))
    }
    if (hasArg(description)) {
        args <- c(args, list(description = quote(description)))
    }
    do.call("flowFrame", args, envir = environment())
}

## used in sort_files_by_component()
sort_component_helper <- function (splits, indices, orders) {
    n <- length(indices)
    if (n == 0) {
        return(indices)
    }
    empty_p <- lapply(indices, function (i) length(splits[[i]]) == 0) %>% unlist
    if (all(empty_p)) {
        return(indices)
    }
    nonempty_inds <- indices[!empty_p]
    next_orders <- orders[-1]
    next_splits <- lapply(splits, function (cur) cur[-1])
    cur_splits <- lapply(indices, function (i) splits[[i]][1]) %>% unlist
    recognized <- if (length(orders) > 0) { orders[[1]] } else { character() }
    lvls <- cur_splits[!(cur_splits %in% recognized)] %>% sort %>% unique %>% {
        as.character(c(recognized, .))
    }
    nonempty_reduced <- Reduce(
        x = lvls,
        init = list(result = list(),
                    remaining = nonempty_inds),
        f = function (cur, level) {
            if (length(cur$remaining) == 0) {
                return(cur)
            }
            matching <- lapply(cur$remaining, function (i) {
                splits[[i]][1] == level
            }) %>% unlist
            matched_sorted <- sort_component_helper(
                next_splits, cur$remaining[matching], next_orders)
            list(result = c(cur$result, matched_sorted),
                 remaining = cur$remaining[!matching])
        })
    stopifnot(length(nonempty_reduced$remaining) == 0)
    c(indices[empty_p], nonempty_reduced$result) %>% unlist
}

sort_files_by_component <- function (strs, split_by, orders = list(),
                                     split_fixed = T, value = T) {
    splits <- if (split_fixed) {
                  strsplit(strs, split_by, fixed = T)
              } else {
                  strsplit(strs, split_by, perl = T)
              }
    indices <- sort_component_helper(splits, 1:length(splits), orders)
    if (value) {
        strs[indices]
    } else {
        indices
    }
}



### Clean fcs data.

numeric_pat <- "^[[:digit:]]+$"
sne_pat <- "[sS][nN][eE]"

is_double_vec <- function (vec) {
    is.vector(vec, mode = "double")
}

is_all_integer <- function (vec) {
    ## this benchmarked as the fastest method
    all(vec == as.integer(vec))
}

## TODO: source this transform and why it's useful
asinh_transform <- function (x) {
    asinh(x / 5)
}

## TODO: make config file to canonicalize column names!
## Return data frame which contains only marker data.
## sample_n: equal sampling
## sample_frac: proportional sampling
fcs_data_cols <- function (
    frame,
    excl_pats = list(numeric_pat, sne_pat),
    excl_preds = list(is_all_integer)
    ) {
    frame %>%
        ## filter out column names by regular expression
        Reduce(init = ., x = excl_pats, f = function (df, pat) {
            df %>% select_if({ !grepl(pat, colnames(.), perl = T) })
        }) %>%
        ## filter out columns by their value or type
        Reduce(init = ., x = excl_preds, f = function(df, pred) {
            ## pred takes a vector and outputs a single logical value
            ## apply pred to each column of the data frame and remove columns
            ## for which it returns TRUE
            df %>% select_if({ !apply(., 2, pred) })
        })
}

shared_markers <- function (frames) {
    ## TODO: find less common columns and check if they're mistakes
    ## TODO: if columns are close but not the same (e.g. levenshtein), show a
    ## warning
    frames %>% lapply(colnames) %>% Reduce(f = intersect)
}



### Manipulate fcs data.

do_tsne <- function (infiles, n,
                     markers = NULL, transform_with = asinh_transform,
                     verbose = T, ...) {
    num_f <- length(infiles)
    suffix <- sprintf("_visne_events_%s.fcs", n)
    outfiles <- gsub("\\.fcs$", suffix, infiles, perl = T)
    if (verbose) {
        cat(sprintf("sampling %s rows per file for %s files...\n",
                    n, num_f))
    }
    sampled_frames <- infiles %>% lapply(function (file) {
        if (verbose) {
            cat(sprintf("reading %s...\n", file))
        }
        read_file(file) %>% sample_n(size = n)
    })
    on_markers <-
        if (!is.null(markers)) {
            markers
        } else {
            if (verbose) {
                cat("guessing markers to join on...\n")
            }
            sampled_frames %>% lapply(fcs_data_cols) %>% shared_markers
        }
    if (verbose) {
        cat(sprintf("joining on markers:\n[%s]\n",
                    paste0(on_markers, collapse = ", ")))
        cat("performing visne...\n")
    }
    ## join into one big table to perform tsne
    joined <- sampled_frames %>% lapply(function (df) {
        ## select the shared columns and squash the raw data
        ## throws error if a frame doesn't have all the markers specified
        df[,on_markers] %>% mutate_all(transform_with)
    }) %>% Reduce(f = rbind)
    tsne_frame <- Rtsne(
        joined, dims = 2, check_duplicates = F, verbose = verbose, ...)
    if (verbose) {
        cat("writing visne results to files...\n")
    }
    ## write to fcs files and check that they exist
    lapply(1:num_f, function (i) {
        sampled_frame <- sampled_frames[[i]]
        outfile <- outfiles[i]
        if (verbose) {
            cat(sprintf("writing %s (%s/%s)...\n", outfile, i, num_f))
        }
        ## find the tsne axes that belong to us
        tsne_start <- ((i - 1) * n) + 1
        tsne_range <- tsne_start:(tsne_start + n - 1)
        ## add tsne axes and convert to flowFrame
        sampled_frame %>%
            mutate(tSNE1 = tsne_frame$Y[tsne_range,1],
                   tSNE2 = tsne_frame$Y[tsne_range,2]) %>%
            make_flowFrame %>% write.FCS(filename = outfile)
    })
    stopifnot(all(file.exists(outfiles)))
    outfiles
}



### Analyze hierarchies of populations in a dataset.

emd_frames <- function (frames, max_rows = 50, max_iterations = 10) {
    mats <- lapply(frames, function (df) {
        df %>% select(c(tSNE1, tSNE2)) %>% slice(1:max_rows) %>%
            as.matrix
    })
    t_emdw <- system.time(ret_emdw <- emdw(mats[[1]], rep(1, max_rows),
                                           mats[[2]], rep(1, max_rows),
                                           max.iter = max_iterations))
    ## t_simplex <- system.time(ret_simplex <- {
    ##     ftm <- t(mats[[1]])
    ##     gtm <- t(mats[[2]])
    ##     n <- max_rows^2
    ##     d_r <- vector(mode = "double", length = n)
    ##     for (i in 1:max_rows) {
    ##         col <- ftm[,i]
    ##         d_r[(1:max_rows) + (max_rows*(i - 1))] <-
    ##             ((gtm - col)^2) %>% colSums %>% sqrt
    ##     }
    ##     lts <- matrix(rep(0.0, max_rows * n), ncol = n)
    ##     gts <- matrix(rep(0.0, max_rows * n), ncol = n)
    ##     for (i in 1:max_rows) {
    ##         lts[i,((1:max_rows) + (max_rows*(i - 1)))] <-
    ##             rep(1.0, max_rows)
    ##         gts[i,(((1:max_rows - 1)*max_rows) + i)] <-
    ##             rep(1.0, max_rows)
    ##     }
    ##     simplex(d_r,
    ##             A1 = lts, b1 = rep(1.0, max_rows),
    ##             A2 = gts, b2 = rep(1.0, max_rows),
    ##             n.iter = max_iterations)
    ## })
    list(t_emdw = t_emdw,
         ret_emdw = ret_emdw
        ## ,
        ##  t_simplex = t_simplex,
        ##  ret_simplex = ret_simplex
         )
}

## TODO: parallelize this!
emd_fcs <- function (files, outfile,
                     max_iterations = 10,
                     verbose = T) {
    stopifnot(!any(duplicated(files)))
    n <- length(files)
    if (verbose) {
        cat(sprintf("reading in %s files...\n", n))
    }
    mats <- lapply(files, function (file) {
        read_file(file) %>% select(c(tSNE1, tSNE2)) %>%
            slice(1:1000) %>% as.matrix
    })
    output <- matrix(vector(mode = "double", length = n * n), nrow = n)
    if (verbose) {
        cat(sprintf("starting pairwise emd on %s files...\n", n))
    }
    for (i in 1:n) {
        if (verbose) {
            cat(sprintf("row %s/%s\n", i, n))
        }
        i_mat <- mats[[i]]
        i_rows <- dim(i_mat)[1]
        if (i > 1) {
            for (j in 1:(i - 1)) {
                if (verbose) {
                    cat(sprintf("column %s/%s\n", j, n))
                }
                output[i,j] <- output[j,i]
            }
        }
        output[i,i] <- 0
        if (verbose) {
            cat(sprintf("column %s/%s\n", i, n))
        }
        if (i < n) {
            for (j in (i + 1):n) {
                if (verbose) {
                    cat(sprintf("column %s/%s\n", j, n))
                }
                j_mat <- mats[[j]]
                j_rows <- dim(j_mat)[1]
                output[i,j] <- emdw(i_mat, rep(1, i_rows),
                                    j_mat, rep(1, j_rows),
                                    max.iter = max_iterations)
            }
        }
    }
    colnames(output) <- files
    rownames(output) <- files
    write.table(output, outfile, sep = ",")
    stopifnot(file.exists(outfile))
    outfile
}

calc_mag_iqr <- function (frame, markers) {
    lapply(markers, function (mark) {
        frame[,mark] %>% { c(median(.), IQR(., type = 2)) }
    }) %>% Reduce(f = rbind) %>%
        set_colnames(c("MAG", "IQR")) %>%
        set_rownames(markers) %>% as.data.frame
}

calc_mem <- function (pop, ref, markers) {
    stopifnot(all(rownames(pop) == markers) &&
              all(rownames(ref) == markers))
    flip_mems <- (pop$MAG - ref$MAG) < 0
    mems <- (abs(pop$MAG - ref$MAG) + (ref$IQR / pop$IQR) - 1) *
        ((-2 * flip_mems) + 1)
    setNames(mems, markers)
}

mem_fcs <- function (files, outfile,
                     markers = NULL,
                     transform_with = asinh_transform,
                     verbose = T) {
    n <- length(files)
    if (verbose) {
        cat(sprintf("performing pairwise MEM on %s files...\n", n))
    }
    frames <- files %>% lapply(function (file) {
        if (verbose) {
            cat(sprintf("reading %s...\n", file))
        }
        read_file(file)
    })
    on_markers <-
        if (!is.null(markers)) {
            markers
        } else {
            if (verbose) {
                cat("guessing markers to join on...\n")
            }
            frames %>% lapply(fcs_data_cols) %>% shared_markers
        }
    if (verbose) {
        cat(sprintf("joining on markers:\n[%s]\n",
                    paste0(on_markers, collapse = ", ")))
        cat("performing MEM...\n")
    }
    marked_pops <- lapply(frames, function (df) {
        df[,on_markers] %>% mutate_all(transform_with)
    })
    stat_dfs <- lapply(marked_pops, function (df) calc_mag_iqr(df, on_markers))
    global_ref <- Reduce(x = marked_pops, f = rbind) %>%
        calc_mag_iqr(on_markers)
    mem_vectors <- lapply(stat_dfs, function (st_df) {
        calc_mem(st_df, global_ref, on_markers)
    })
    output <- matrix(vector("double", length = n * n), nrow = n)
    if (verbose) {
        cat(sprintf("starting pairwise MEM RMSD on %s files...\n", n))
    }
    for (i in 1:n) {
        if (verbose) {
            cat(sprintf("row %s/%s\n", i, n))
        }
        i_vt <- mem_vectors[[i]]
        for (j in 1:n) {
            if (verbose) {
                cat(sprintf("col %s/%s\n", j, n))
            }
            j_vt <- mem_vectors[[j]]
            output[i,j] <- (i_vt - j_vt) ^ 2 %>% sum %>% sqrt
        }
    }
    colnames(output) <- files
    rownames(output) <- files
    write.table(output, outfile, sep = ",")
    stopifnot(file.exists(outfile))
    outfile
}

## look at xtabs/ftable/table() and summary/summarize/aggregate/group_by()



### Pull data from cytobank.

## TRUE unless file exists, is readable, and has the right size and contents
invalid_fcs_dl <- function (fcs_info) {
    fcs_info %$%
        {
            (file.exists(filename) &
             file.access(filename, mode = 4) == 0 &
             fileSize == file.size(filename) &
             md5sum == md5sum(filename))
        } %>%
        ## any NA or F means file is inaccessible or corrupted
        ## (e.g. by an interrupted download)
        { is.na(.) | !(.) }
}

## return fcs filenames
download_all_fcs <- function (session, exp_id, verbose = T) {
    fcs_info <- fcs_files.list(session, exp_id) %>%
        ## fcs_files.list returns a data frame with list columns; undo that
        mutate_all(unlist) %>%
        select(id, filename, md5sum, fileSize) %>%
        mutate(humansize = humanReadable(fileSize, width = 4))
    to_dl <- fcs_info[invalid_fcs_dl(fcs_info),]
    size_all <- sum(to_dl$fileSize)
    size_cur <- 0
    if (verbose) {
        cat(sprintf("total (unzipped) size of files to download: %s\n",
                    humanReadable(size_all, width = 4)))
    }
    num_to_dl <- dim(to_dl)[1]
    if (num_to_dl != 0) {
        for (i in 1:num_to_dl) {
            row <- to_dl[i,]
            if (verbose) {
                cat(sprintf("%s downloaded (%s%% done)\n",
                            humanReadable(size_cur, width = 4),
                            round(size_cur / size_all * 100, digits = 1)))
                cat(sprintf("downloading %s (unzipped size: %s, id: %s)...\n",
                            row$filename, row$humansize, row$id))
            }
            taken <- system.time(
                fcs_files.download_zip(session, exp_id, row$id) %>% unzip)
            cat(sprintf("download of %s took %s seconds\n",
                        row$filename, round(taken[3], 2)))
            if (any(invalid_fcs_dl(row))) {
                stop(sprintf(paste("download of %s is corrupt",
                                   "-- check your internet connection"),
                             row$filename))
            }
            size_cur <- size_cur + row$fileSize
        }
    }
    fcs_info$filename
}

download_fcs_cytobank <- function (experiments, verbose = T) {
    if (length(experiments) == 0) { return(character()) }
    library(getPass)
    session <- authenticate(readline("site (<site>.cytobank.org): "),
                            readline("username: "),
                            getPass(forcemask = T))
    Reduce(f = c, init = character(),
           x = lapply(experiments, function (exp_id) {
               download_all_fcs(session, exp_id, verbose = verbose)
           }))
}

download_gates <- function (session, exp_id) {
    gates.gatingML_download(session, exp_id)
}

apply_gates <- function (gates_xml, fcs_files) {
    cytobank2GatingSet(gates_xml, fcs_files)
}



### Parsing cytobank output.
.xpath <- function (doc, xpath_str = ".", node = doc, fun = NULL) {
    XML::xpathSApply(doc = node, path = xpath_str, fun = fun,
                     namespaces = xmlNamespaceDefinitions(doc, simplify = T))
}

## doc <- xmlParse("./allie-paper/CytExp_22899_Gates_v1.xml")
## fcs <- list.files(path = "./allie-paper/", pattern = "fcs$", full.names = T)
## gates <- .parse_rectangle_gate(doc)
## processed_fcs <- lapply(fcs, function (file) {
##     Reduce(init = read_file(file), x = gates, f = function (acc, cur) {
##         add_to(cur, acc)
##     })
## })
.parse_rectangle_gate <- function(xml) {
    .xpath(xml, "/gating:Gating-ML/gating:RectangleGate") %>%
        lapply(function (rect_gate_xml) {
            new("RectangleGate", rect_gate_xml, xml)
        })
}

.parse_polygon_gate <- function (gate_node) {
    gate_node
}

.data_transformation_fun_dict <- list(
    Tr_Arcsinh_5 = asinh_transform
)

.anon <- function (block, env = parent.frame()) {
    sb <- substitute(block)
    eval(bquote(function (.) { eval(.(sb), .(env)) }))
}

.explode <- function (result, pred = is.null) {
    stopifnot(isTRUE(!((match.fun(pred))(result))))
    result
}

.get_single <- function (lst) {
    .explode(lst, function (l) length(l) != 1) %>% .[[1]]
}

setGeneric('add_to', function (gate, fcs) stop(paste(gate, fcs)))

setClass('Gate',
         slots = c(gate_name = 'character', gate_id = 'character'))
setMethod('initialize', 'Gate', function(.Object, gate_xml, doc) {
    .Object <- callNextMethod(.Object)
    gate_name <- .xpath(doc, 'data-type:custom_info/cytobank/name/text()',
                        gate_xml, xmlValue) %>% .get_single
    gate_id <- xmlGetAttr(gate_xml, 'gating:id')
    .Object@gate_name <- gate_name
    .Object@gate_id <- gate_id
    .Object
})

setClass('RectangleGate',
         slots = c(constraints = 'list'),
         contains = 'Gate')
setMethod(
    'initialize', 'RectangleGate',
    function (.Object, rect_gate_xml, doc) {
        .Object <- callNextMethod(.Object, rect_gate_xml, doc)
        constraints <- .xpath(doc, 'gating:dimension', rect_gate_xml) %>%
            lapply(function (dim_xml) {
                list(transform_name = xmlGetAttr(dim_xml,
                                                 'gating:transformation-ref'),
                    marker = .xpath(
                        doc, 'data-type:fcs-dimension/@data-type:name',
                        dim_xml),
                    min = xmlGetAttr(dim_xml, 'gating:min'),
                    max = xmlGetAttr(dim_xml, 'gating:max'))
            })
        .Object@constraints <- constraints
        .Object
    })
setMethod('add_to', c(gate = 'RectangleGate', fcs = 'data.frame'),
          function (gate, fcs) {
              fcs[,gate@gate_name] <- TRUE
              Reduce(x = gate@constraints, init = fcs, f = function (acc, cur) {
                  tr_f <- match.fun(
                      .data_transformation_fun_dict[[cur$transform_name]])
                  acc[,gate@gate_name] <- acc[,gate@gate_name] &
                      acc[,cur$marker] %>% tr_f %>% {
                          (. >= cur$min) & (. <= cur$max)
                      }
                  acc
              })
          })

## setClassUnion(
##     "Gate", c("RectangleGate", "PolygonGate", "BooleanGate", "QuadrantGate"))

.parse_quadrant_gate <- function (gatingML, files) {
    doc <- xmlParse(gatingML)
    quadrant_gate_nodes <- .xpath(doc, "/gating:Gating-ML/gating:QuadrantGate")
    lapply(quadrant_gate_nodes, function (node) {
        ## name <-
    })
}

.parse_boolean_gate <- function (gate_node) {
    gate_node
}

.gate_parse_dispatch <- list(
    ## TODO: check if we're missing any gate types in the parse function!
    PolygonGate = .parse_polygon_gate,
    RectangleGate = .parse_rectangle_gate,
    QuadrantGate = .parse_quadrant_gate,
    BooleanGate = .parse_boolean_gate
)

parse_gatingml <- function (gatingMLFile) {
    doc <- xmlParse(gatingMLFile)
    lapply(.gate_parse_dispatch, function (p) p(doc)) %>% Reduce(f = c)
}

## .gate_apply_dispatch <- list(
##     PolygonGate = .apply_polygon_gate,
##     RectangleGate = .apply_rectangle_gate,
##     QuadrantGate = .apply_quadrant_gate,
##     BooleanGate = .apply_boolean_gate
## )

apply_parsed_gates <- function (parsed_gates, cyto_files) {
    lapply(cyto_files, read_file) %>%
        Reduce(x = parsed_gates, init = ., f = function (frame_list, gate) {
            lapply(frame_list, function (cyto) apply_gate(gate, cyto))
        })
}

parse_gates_pops <- function (xml) {
    doc <- xmlParse(xml)
    doc_ns <- xmlNamespaceDefinitions(doc, simplify = T)
    all_gate_types <- xpathSApply(
        doc, "/gating:Gating-ML/gating:*[contains(., 'Gate')]", xmlName) %>%
        unique
    recognized_gate_types <- names(gate_analysis_dispatch)
    stopifnot(all(!duplicated(recognized_gate_types)) &&
              setequal(all_gate_types, recognized_gate_types))

}



### Comparative clustering.

daniel_read_fcs <- function (fname) {
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

f_measure_clusters <- function (frame) {
    cls <- frame$cluster %>% unique %>% sort
    k <- length(cls)
    stopifnot(all(cls == 1:k))
    lapply(cls, function (cl_ind) {
        mode_pop <- frame[frame$cluster == cl_ind,] %>% .$manual_gates %>%
            table %>% sort(decreasing = T) %>%
            names %>% as.integer %>% .[1]
        n_relevant <- sum(frame$manual_gates == mode_pop &
                          frame$cluster == cl_ind)
        precision <- n_relevant / sum(frame$cluster == cl_ind)
        recall <- n_relevant / sum(frame$manual_gates == mode_pop)
        2 * precision * recall / (precision + recall)
    }) %>% unlist %>% mean
}

cluster_kmeans <- function (frame, k) {
    clustering <- frame %>% select(c(tSNE1, tSNE2)) %>% kmeans(centers = k)
    frame %>% mutate(cluster = clustering$clust)
}

cluster_spade <- function (files) {
    k <- lapply(files, read_file) %>% Reduce(f = bind_rows) %>% num_term_pops
    cat(sprintf("k = %s\n", k))
    suppressWarnings(
        SPADE.driver(files = files, cluster_cols = c("tSNE1", "tSNE2"),
                     transforms = NULL, k = k,
                     downsampling_target_percent = .01))
    clustered <- sprintf("%s.density.fcs.cluster.fcs", files)
    lapply(clustered, read_fcs) %>% Reduce(f = bind_rows)
}

##
##
##

## ### Run viSNE on a set of files, compute pairwise Earth Mover's Distance (EMD),
## ### and generate a heatmap.
## ## Written by Danny McClanahan, Irish Lab June 2017.
## ## <danieldmcclanahan@gmail.com>

## ### Configure
## ## `color_palette`: Color palette for heatmap.
## color_palette <- colorRampPalette(
##     c("#24658C", "#5CBAA7", "#9ED2A4", "#E2E998", "#FBF7BF", "#FDDC86",
##       "#F8A05A", "#EF6342", "#D43E4F")
## )(n = 50)
## ## `emd_outfile`: filename for output csv with EMD pairwise comparisons
## emd_outfile <- "emd_out.csv"
## ## `emd_heatmap_outfile`: output file containing EMD heatmap as a pdf
## emd_heatmap_outfile <- "heatmap_emd.pdf"
## ## `mem_outfile`: filename for output csv with MEM RMSD pairwise comparisons
## mem_outfile <- "mem_out.csv"
## ## `mem_heatmap_outfile`: output file containing MEM heatmap as a pdf
## mem_heatmap_outfile <- "heatmap_mem.pdf"
## ## Set working directory, if you need to.
## setwd(".")

## ## data_files: Char vector of files to read data from.
## ##   Files can be binary fcs files or text files with headers. ".txt" files will
## ##   be read as TSV (sep = "\t"), and ".csv" files as CSV (sep = ",").
## ##
## ##   NOTE: Files should have two tSNE axes labeled tSNE1 and tSNE2!!!
## ##
## ##   If the file is a text file (.txt or .csv) and starts with a blank line,
## ##   this script will recognize that and skip the initial blank line.
## ##
## ##   You can manually set data_files as well. For example:
## ##   > data_files <- c("file1.fcs", "file2.fcs")
## data_files <- list.files(pattern = "\\.fcs$", ignore.case = TRUE,
##                          all.files = TRUE, full.names = TRUE, recursive = FALSE,
##                          no.. = TRUE)

## ## sort_files_by_component(): splits filenames into pieces and sorts them.
## ##   sort_files_by_component() uses the string in `split_by` to break filenames
## ##   into components.
## ##
## ##   Splitting "A:B:C" by ":" returns c("A", "B", "C")). We say that the string
## ##   "A:B:C" has components "A" at index 1, "B" at index 2, and "C" at 3.
## ##
## ##   `orders` is a list of character vectors. If a component shows up here at
## ##   the correct index, its file is ordered before other files. Among files
## ##   which have matching components at each index, the order is determined by
## ##   the order in the character vector at that index.
## ##
## ##   Example:
## ##
## ##   > data_files <- c("MB004_6m_panel2.fcs", "MB004_3wk_panel2.fcs",
## ##                     "MB004_12wk_panel2.fcs", "MB004_pre_panel2.fcs",
## ##                     "MB005_12wk_panel2.fcs")
## ##   > sort_files_by_component(
## ##       data_files,
## ##       split_by = "_",
## ##       orders = list(c("MB005"), c("pre", "3wk", "12wk", "6m")))
## ##   [1] "MB005_12wk_panel2.fcs" "MB004_pre_panel2.fcs"  "MB004_3wk_panel2.fcs"
## ##   [4] "MB004_12wk_panel2.fcs" "MB004_6m_panel2.fcs"
## ##
## ##   "MB005" sorts before "MB004" in the first component above because "MB005"
## ##   is mentioned in `orders` at index 1, but "MB004" isn't. If we did it again
## ##   with `orders` = list(c(), c("pre", "3wk", "12wk", "6m")), the first
## ##   component would be sorted alphabetically, so "MB005" would come after
## ##   "MB004".
## ##
## ##   With `split_by` = "" and `orders` = list(), `data_files` is simply sorted
## ##   alphabetically by file name.
## data_files_sorted <- sort_files_by_component(
##     data_files,
##     split_by = "_",
##     orders = list(c(), c("pre", "3wk", "12wk", "6m")))

## ## emd_fcs(): Run pairwise EMD on input files and produce CSV.
## ##   The resuts are stored in `emd_outfile`.
## ##
## ##   `max_iterations` is the number of iterations to perform when computing EMD.
## ##   Increasing this value *typically* does not change the result at all.
## emd_fcs(data_files_sorted, emd_outfile, max_iterations = 10)
## ## emd_outfile has row names in column 1
## emd_matrix <- as.matrix(read.csv(emd_outfile, row.names = 1))

## pdf(emd_heatmap_outfile)
## heatmap.2(emd_matrix, Rowv = F, Colv = F, dendrogram = "none",
##           col = color_palette, trace = "none", density.info = "none")
## dev.off()

## ## mem_fcs(): Run pairwise MEM RMSD on input files and produce CSV.
## ##   The resuts are stored in `mem_outfile`.
## mem_fcs(data_files_sorted, mem_outfile)
## ## mem_outfile has row names in column 1
## mem_matrix <- as.matrix(read.csv(mem_outfile, row.names = 1))

## pdf(mem_heatmap_outfile)
## heatmap.2(
##     mem_matrix, col = color_palette,
##     Rowv = F, Colv = F, dendrogram = "none", trace = "none",
##     density.info = "none")
## dev.off()
