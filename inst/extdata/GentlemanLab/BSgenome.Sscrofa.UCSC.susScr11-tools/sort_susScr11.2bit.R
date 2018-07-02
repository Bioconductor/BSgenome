library(Biostrings)
library(rtracklayer)
susScr11 <- import("susScr11.2bit")

seqlevels1 <- paste0("chr", c(1:18, "X", "Y", "M"))

tmp <- CharacterList(strsplit(names(susScr11), "_"))
npart <- elementNROWS(tmp)
stopifnot(all(npart <= 4L))

idx1 <- which(npart == 1L)
stopifnot(setequal(names(susScr11)[idx1], seqlevels1))

idx2 <- which(npart == 2L)
stopifnot(length(idx2) == 0L)

idx4 <- which(npart == 4L)
m4 <- matrix(unlist(tmp[idx4]), ncol=4L, byrow=TRUE)
stopifnot(all(m4[ , 2L] == "NW"))
stopifnot(all(m4[ , 4L] == "random"))
chrom <- match(m4[ , 1L], seqlevels1)
stopifnot(!anyNA(chrom))
oo4 <- order(chrom, m4[ , 3L])
seqlevels4 <- names(susScr11)[idx4[oo4]]

idx3 <- which(npart == 3L)
m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
stopifnot(all(m3[ , 1L] == "chrUn"))
stopifnot(all(m3[ , 2L] == "NW"))
oo3 <- order(m3[ , 3L])
seqlevels3 <- names(susScr11)[idx3[oo3]]

susScr11_seqlevels <- c(seqlevels1, seqlevels4, seqlevels3)

export(susScr11[susScr11_seqlevels], "susScr11.sorted.2bit")

