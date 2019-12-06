library(Biostrings)
library(rtracklayer)
calJac3 <- import("calJac3.2bit")

seqlevels1 <- paste0("chr", c(1:22, "X", "Y"))

tmp <- CharacterList(strsplit(names(calJac3), "_"))
npart <- lengths(tmp)
stopifnot(all(npart <= 3L))

idx1 <- which(npart == 1L)
stopifnot(setequal(names(calJac3)[idx1], seqlevels1))

idx3 <- which(npart == 3L)
m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
m31 <- match(m3[ , 1L], seqlevels1)
stopifnot(!anyNA(m31))
stopifnot(all(m3[ , 3L] == "random"))
oo3 <- order(m31, m3[ , 2L])
seqlevels3 <- names(calJac3)[idx3[oo3]]

idx2 <- which(npart == 2L)
m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
stopifnot(all(m2[ , 1L] == "chrUn"))
oo2 <- order(m2[ , 2L])
seqlevels2 <- names(calJac3)[idx2[oo2]]

calJac3_seqlevels <- c(seqlevels1, seqlevels3, seqlevels2)

export(calJac3[calJac3_seqlevels], "calJac3.sorted.2bit")

