### Routines for manipulating CyToF data and interacting with Cytobank.
## Written by Danny McClanahan, Irish Lab June 2017.
## <danieldmcclanahan@gmail.com>



### Libraries
### Keep the libraries up to date
library(flowCore, warn.conflicts = F)
library(tools)
library(digest)
library(boot)
library(XML)
library(gplots)
library(emdist, warn.conflicts = F)
library(spade, quietly = T, warn.conflicts = F, verbose = F)
library(magrittr, warn.conflicts = F)
library(dplyr, warn.conflicts = F)



### Read/write different representations of flow data.

## NOTE: we only grab a single dataset from each input file, even if more
## exist. Multiple datasets in an fcs file is deprecated and "implementors are
## being discouraged to do so." it may exist, but shouldn't be supported.

.read_fcs_flowFrame <- function (fname) {
    read.FCS(fname, transformation = NULL, truncate_max_range = F)
}

.read_fcs_cyto_frame <- function (fname) {
    .read_fcs_flowFrame(fname) %>% exprs %>% as.data.frame
}

.read_text_cyto_frame <- function (fname, ...) {
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
        fcs = .read_fcs_cyto_frame(fname),
        csv = .read_text_cyto_frame(fname, sep = ","),
        txt = .read_text_cyto_frame(fname, sep = "\t"),
        stop(sprintf("unrecognized extension '%s' for file '%s'",
                     ext, fname)))
}

## TODO: check validity of data in df? against params/description?
.make_flowFrame <- function (exprs, parameters, description) {
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
.sort_component_helper <- function (splits, indices, orders) {
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
            matched_sorted <- .sort_component_helper(
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
    indices <- .sort_component_helper(splits, 1:length(splits), orders)
    if (value) {
        strs[indices]
    } else {
        indices
    }
}



### Clean fcs data.

.numeric_pat <- "^[[:digit:]]+$"
.sne_pat <- "[sS][nN][eE]"

.is_double_vec <- function (vec) {
    is.vector(vec, mode = "double")
}

.is_all_integer <- function (vec) {
    ## this benchmarked as the fastest method
    all(vec == as.integer(vec))
}

## TODO: get a citation for this transform and justify it. get comparisons
asinh_transform <- function (x) { asinh(x / 5) }

## TODO: make config file to canonicalize column names!
## Return data frame which contains only marker data.
.fcs_data_cols <- function (
    frame,
    excl_pats = list(.numeric_pat, .sne_pat),
    excl_preds = list(.is_all_integer)
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


### Analyze hierarchies of populations in a dataset.

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

.calc_mag_iqr <- function (frame, markers) {
    lapply(markers, function (mark) {
        frame[,mark] %>% { c(median(.), IQR(., type = 2)) }
    }) %>% Reduce(f = rbind) %>%
        set_colnames(c("MAG", "IQR")) %>%
        set_rownames(markers) %>% as.data.frame
}

.calc_mem <- function (pop, ref, markers) {
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
            frames %>% lapply(.fcs_data_cols) %>% shared_markers
        }
    if (verbose) {
        cat(sprintf("joining on markers:\n[%s]\n",
                    paste0(on_markers, collapse = ", ")))
        cat("performing MEM...\n")
    }
    marked_pops <- lapply(frames, function (df) {
        df[,on_markers] %>% mutate_all(transform_with)
    })
    stat_dfs <- lapply(marked_pops, function (df) .calc_mag_iqr(df, on_markers))
    global_ref <- Reduce(x = marked_pops, f = rbind) %>%
        .calc_mag_iqr(on_markers)
    mem_vectors <- lapply(stat_dfs, function (st_df) {
        .calc_mem(st_df, global_ref, on_markers)
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



### Wrapper functions for xpath selectors.
.xpath <- function (doc, xpath_str = ".", node = doc, fun = NULL) {
    doc_ns <- xmlNamespaceDefinitions(doc, simplify = T)
    XML::xpathSApply(
        doc = node, path = xpath_str, fun = fun, namespaces = doc_ns)
}



### Parsing cytobank gatingML xml files.

## .gen_generic <- function (name, ) {
## }

setGeneric("add_gate", function (gate, fcs) {

})

setClass("Gate",
         slots = c(gate_name = "character", gate_id = "character"))
setMethod("initialize", "Gate", function(.Object, gate_xml, doc) {
    .Object <- callNextMethod(.Object)
    gate_name <- .xpath(doc, "data-type:custom_info/cytobank/name/text()",
                        gate_xml, xmlValue) %>% .get_single
    gate_id <- xmlGetAttr(gate_xml, "gating:id")
    .Object@gate_name <- gate_name
    .Object@gate_id <- gate_id
    .Object
})

.data_transformation_fun_dict <- list(
    Tr_Arcsinh_5 = asinh_transform
)

setClass("RectangleGate",
         slots = c(constraints = "list"),
         contains = "Gate")
setMethod(
    "initialize", "RectangleGate",
    function (.Object, rect_gate_xml, doc) {
        .Object <- callNextMethod(.Object, rect_gate_xml, doc)
        constraints <- .xpath(doc, "gating:dimension", rect_gate_xml) %>%
            lapply(function (dim_xml) { list(
                transform_name = .xpath(
                    doc, "@gating:transformation-ref", dim_xml),
                marker = .xpath(
                    doc, "data-type:fcs-dimension/@data-type:name", dim_xml),
                min = .xpath(doc, "@gating:min", dim_xml),
                max = .xpath(doc, "gating:max", dim_xml))
            })
        .Object@constraints <- constraints
        .Object
    })
setMethod(
    "add_gate", c(gate = "RectangleGate", fcs = "data.frame"),
    function (gate, fcs) {
        fcs[,gate@gate_name] <- TRUE
        Reduce(
            init = fcs,
            x = gate@constraints,
            f = function (acc, constr) {
                tr_f <- match.fun(
                    .data_transformation_fun_dict[[constr$transform_name]])
                satisfies_constraint <- acc[,constr$marker] %>% tr_f %>% {
                    (. >= constr$min) & (. <= constr$max)
                }
                acc[,gate@gate_name] <-
                    acc[,gate@gate_name] & satisfies_constraint
                acc
            })
    })

## setClassUnion(
##     "Gate", c("RectangleGate", "PolygonGate", "BooleanGate", "QuadrantGate"))

## doc <- xmlParse("./allie-paper/CytExp_22899_Gates_v1.xml")
## fcs <- list.files(path = "./allie-paper/", pattern = "fcs$", full.names = T)
## gates <- .parse_rectangle_gates(doc)
## processed_fcs <- lapply(fcs, function (file) {
##     Reduce(init = read_file(file), x = gates, f = function (acc, cur) {
##         add_gate(cur, acc)
##     })
## })

.parse_rectangle_gates <- function (xml) {
    ## "data-type:custom_info/cytobank/fcs_file_filename/text()"
    .xpath(xml, "/gating:Gating-ML/gating:RectangleGate") %>%
        lapply(function (rect_gate_node) {
            new("RectangleGate", rect_gate_node, xml)
        })
}

.parse_polygon_gates <- function (xml) {
    xml
}

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

.parse_quadrant_gates <- function (xml) {
    quadrant_gate_nodes <- .xpath(xml, "/gating:Gating-ML/gating:QuadrantGate")
    lapply(quadrant_gate_nodes, function (quadrant_gate_node) {
    })
}

.parse_boolean_gates <- function (xml) {
    xml
}

.gate_parse_dispatch <- list(
    ## TODO: check if there are node names matching /^gating:.*Gate$/ that
    ## aren't in this list
    PolygonGate = .parse_polygon_gates,
    RectangleGate = .parse_rectangle_gates,
    QuadrantGate = .parse_quadrant_gates,
    BooleanGate = .parse_boolean_gates
)
