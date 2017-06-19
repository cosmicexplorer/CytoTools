library(flowCore)
library(dplyr, warn.conflicts = F)

read_clean_fcs <- function(fname, clean = T, change_names = F) {
    flow_fcs <- read.FCS(fname, alter.names = change_names)
    flow_frame <- as.data.frame(exprs(flow_fcs))

    if (!clean) { flow_frame }
    else { flow_frame %>% na.omit %>% distinct(.keep_all = T) }
}

squish_expression <- Vectorize(function(x) asinh(x / 5))
inv_squish_expression <- Vectorize(function(y) 5 * sinh(y))

add_squished <- function (frame, channels,
                          sep = "S-", squish = squish_expression) {
    squished_names <- paste(sep, colnames(frame)[which(channels)], sep = "")
    squished <- squish(frame[,channels])
    colnames(squished) <- squished_names
    cbind(frame[,!channels], frame[,channels], squished)
}
