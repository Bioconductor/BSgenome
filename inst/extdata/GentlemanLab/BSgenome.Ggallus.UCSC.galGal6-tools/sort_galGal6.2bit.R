library(Biostrings)
library(rtracklayer)
galGal6 <- import("galGal6.2bit")

seqlevels1 <- paste0("chr", c(1:28, 30:33, "M", "W", "Z"))

tmp <- CharacterList(strsplit(names(galGal6), "_"))
npart <- lengths(tmp)
stopifnot(all(npart %in% c(1L, 3L, 4L)))

idx1 <- which(npart == 1L)
stopifnot(setequal(names(galGal6)[idx1], seqlevels1))

idx4 <- which(npart == 4L)
m4 <- matrix(unlist(tmp[idx4]), ncol=4L, byrow=TRUE)
m41 <- match(m4[ , 1L], seqlevels1)
stopifnot(!anyNA(m41))
stopifnot(all(m4[ , 2L] == "NW"))
stopifnot(all(m4[ , 4L] == "random"))
oo4 <- order(m41, m4[ , 3L])
seqlevels4 <- names(galGal6)[idx4[oo4]]

idx3 <- which(npart == 3L)
m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
stopifnot(all(m3[ , 1L] == "chrUn"))
stopifnot(all(m3[ , 2L] == "NW"))
oo3 <- order(m3[ , 3L])
seqlevels3 <- names(galGal6)[idx3[oo3]]

galGal6_seqlevels <- c(seqlevels1, seqlevels4, seqlevels3)

export(galGal6[galGal6_seqlevels], "galGal6.sorted.2bit")

