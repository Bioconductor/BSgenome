library(Biostrings)
library(rtracklayer)
rheMac8 <- import("rheMac8.2bit")

seqlevels1 <- paste0("chr", c(1:20, "X", "Y", "M"))

tmp <- CharacterList(strsplit(names(rheMac8), "_"))
npart <- elementNROWS(tmp)
stopifnot(all(npart <= 3L))

idx1 <- which(npart == 1L)
stopifnot(setequal(names(rheMac8)[idx1], seqlevels1))

idx2 <- which(npart == 2L)
stopifnot(length(idx2) == 0L)

idx3 <- which(npart == 3L)
m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
stopifnot(all(m3[ , 1L] == "chrUn"))
stopifnot(all(m3[ , 2L] == "NW"))
oo3 <- order(m3[ , 3L])
seqlevels3 <- names(rheMac8)[idx3[oo3]]

rheMac8_seqlevels <- c(seqlevels1, seqlevels3)

export(rheMac8[rheMac8_seqlevels], "rheMac8.sorted.2bit")

