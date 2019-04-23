library(Biostrings)
library(rtracklayer)
monDom5 <- import("monDom5.2bit")

monDom5_seqlevels <- paste0("chr", c(1:8, "X", "M", "Un"))

stopifnot(setequal(names(monDom5), monDom5_seqlevels))

export(monDom5[monDom5_seqlevels], "monDom5.sorted.2bit")

