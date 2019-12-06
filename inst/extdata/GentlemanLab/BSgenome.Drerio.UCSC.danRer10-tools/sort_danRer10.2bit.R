library(Biostrings)
library(rtracklayer)
danRer10 <- import("danRer10.2bit")

seqlevels1 <- paste0("chr", c(1:25, "M"))

tmp <- CharacterList(strsplit(names(danRer10), "_"))
npart <- lengths(tmp)
stopifnot(all(npart <= 2L))

idx1 <- which(npart == 1L)
stopifnot(setequal(names(danRer10)[idx1], seqlevels1))

idx2 <- which(npart == 2L)
m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
stopifnot(all(m2[ , 1L] == "chrUn"))
oo2 <- order(m2[ , 2L])
seqlevels2 <- names(danRer10)[idx2[oo2]]

danRer10_seqlevels <- c(seqlevels1, seqlevels2)

export(danRer10[danRer10_seqlevels], "danRer10.sorted.2bit")

