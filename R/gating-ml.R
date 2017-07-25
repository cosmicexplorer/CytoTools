#' @import magrittr


### Wrapper functions for xpath selectors.

xpath <- function (doc, xpath_str = ".", node = doc, fun = NULL) {
    doc_ns <- XML::xmlNamespaceDefinitions(doc, simplify = TRUE)
    XML::xpathSApply(
        doc = node, path = xpath_str, fun = fun, namespaces = doc_ns)
}



### Parsing cytobank gatingML xml files.

## gen_generic <- function (name, ) {
## }

setGeneric("add_gate", function (gate, cyto_df) {

})

setClass("Gate",
         slots = c(gate_name = "character", gate_id = "character"))
setMethod("initialize", "Gate", function(.Object, gate_xml, doc) {
    .Object <- callNextMethod(.Object)
    gate_name <- xpath(doc, "data-type:custom_info/cytobank/name/text()",
                       gate_xml, xmlValue) %>% get_single
    gate_id <- xmlGetAttr(gate_xml, "gating:id")
    .Object@gate_name <- gate_name
    .Object@gate_id <- gate_id
    .Object
})

data_transformation_fun_dict <- list(
    Tr_Arcsinh_5 = asinh_transform
)

setClass("RectangleGate",
         slots = c(constraints = "list"),
         contains = "Gate")
setMethod(
    "initialize", "RectangleGate",
    function (.Object, rect_gate_xml, doc) {
        .Object <- callNextMethod(.Object, rect_gate_xml, doc)
        constraints <- xpath(doc, "gating:dimension", rect_gate_xml) %>%
            lapply(function (dim_xml) { list(
                transform_name = xpath(
                    doc, "@gating:transformation-ref", dim_xml),
                marker = xpath(
                    doc, "data-type:fcs-dimension/@data-type:name", dim_xml),
                min = xpath(doc, "@gating:min", dim_xml),
                max = xpath(doc, "gating:max", dim_xml))
            })
        .Object@constraints <- constraints
        .Object
    })
setMethod(
    "add_gate", c(gate = "RectangleGate", cyto_df = "data.frame"),
    function (gate, cyto_df) {
        cyto_df[,gate@gate_name] <- TRUE
        Reduce(
            init = cyto_df,
            x = gate@constraints,
            f = function (acc, constr) {
                tr_f <- match.fun(
                    data_transformation_fun_dict[[constr$transform_name]])
                satisfies_constraint <- acc[,constr$marker] %>% tr_f %>% {
                    (. >= constr$min) & (. <= constr$max)
                }
                acc[,gate@gate_name] <-
                    acc[,gate@gate_name] & satisfies_constraint
                acc
            })
    })

## example:
## https://github.com/RGLab/CytoML/blob/trunk/R/read.gatingML.cytobank.R

## setClassUnion(
##     "Gate", c("RectangleGate", "PolygonGate", "BooleanGate", "QuadrantGate"))

## doc <- xmlParse("./allie-paper/CytExp_22899_Gates_v1.xml")
## cyto_df <- list.files(
##     path = "./allie-paper/", pattern = "fcs$", full.names = TRUE)
## gates <- parse_rectangle_gates(doc)
## processed_cyto_frames <- lapply(cyto_df, function (cyto_data_file) {
##     Reduce(init = read_file(cyto_data_file), x = gates,
##            f = function (acc, cur) {
##                add_gate(cur, acc)
##            })
## })

parse_rectangle_gates <- function (xml) {
    ## "data-type:custom_info/cytobank/fcs_file_filename/text()"
    xpath(xml, "/gating:Gating-ML/gating:RectangleGate") %>%
        lapply(function (rect_gate_node) {
            new("RectangleGate", rect_gate_node, xml)
        })
}

parse_polygon_gates <- function (xml) {
    xml
}

parse_quadrant_gates <- function (xml) {
    quadrant_gate_nodes <- xpath(xml, "/gating:Gating-ML/gating:QuadrantGate")
    lapply(quadrant_gate_nodes, function (quadrant_gate_node) {
    })
}

parse_boolean_gates <- function (xml) {
    xml
}

gate_parse_dispatch <- list(
    ## TODO: check if there are node names matching /^gating:.*Gate$/ that
    ## aren't in this list
    PolygonGate = parse_polygon_gates,
    RectangleGate = parse_rectangle_gates,
    QuadrantGate = parse_quadrant_gates,
    BooleanGate = parse_boolean_gates
)
