library(Biostrings)
library(rtracklayer)
bosTau9 <- import("bosTau9.2bit")

seqlevels1 <- paste0("chr", c(1:29, "X", "M"))

tmp <- CharacterList(strsplit(names(bosTau9), "_"))
npart <- lengths(tmp)
stopifnot(all(npart %in% c(1L, 3L)))

idx1 <- which(npart == 1L)
stopifnot(setequal(names(bosTau9)[idx1], seqlevels1))

idx3 <- which(npart == 3L)
m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
stopifnot(all(m3[ , 1L] == "chrUn"))
stopifnot(all(m3[ , 2L] == "NW"))
oo3 <- order(m3[ , 3L])
seqlevels3 <- names(bosTau9)[idx3[oo3]]

bosTau9_seqlevels <- c(seqlevels1, seqlevels3)

export(bosTau9[bosTau9_seqlevels], "bosTau9.sorted.2bit")

