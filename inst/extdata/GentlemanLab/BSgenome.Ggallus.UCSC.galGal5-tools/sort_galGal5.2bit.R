library(Biostrings)
library(rtracklayer)
galGal5 <- import("galGal5.2bit")

seqlevels1 <- paste0("chr", c(1:28, 30:33, "M", "W", "Z", "LGE64"))

tmp <- CharacterList(strsplit(names(galGal5), "_"))
npart <- lengths(tmp)
stopifnot(all(npart %in% c(1L, 3L, 4L)))

idx1 <- which(npart == 1L)
stopifnot(setequal(names(galGal5)[idx1], seqlevels1))

idx4 <- which(npart == 4L)
m4 <- matrix(unlist(tmp[idx4]), ncol=4L, byrow=TRUE)
m41 <- match(m4[ , 1L], seqlevels1)
stopifnot(!anyNA(m41))
stopifnot(all(m4[ , 2L] == "NT"))
stopifnot(all(m4[ , 4L] == "random"))
oo4 <- order(m41, m4[ , 3L])
seqlevels4 <- names(galGal5)[idx4[oo4]]

idx3 <- which(npart == 3L)
m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
stopifnot(all(m3[ , 1L] == "chrUn"))
stopifnot(all(m3[ , 2L] == "NT"))
oo3 <- order(m3[ , 3L])
seqlevels3 <- names(galGal5)[idx3[oo3]]

galGal5_seqlevels <- c(seqlevels1, seqlevels4, seqlevels3)

export(galGal5[galGal5_seqlevels], "galGal5.sorted.2bit")

