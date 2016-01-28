library(Biostrings)
library(rtracklayer)
taeGut2 <- import("taeGut2.2bit")

seqlevels1 <- paste0("chr", c(1, "1A", "1B", 2:4, "4A", 5:28, "LGE22", "LG2", "LG5", "Z", "M"))

tmp <- CharacterList(strsplit(names(taeGut2), "_"))
npart <- elementNROWS(tmp)
stopifnot(all(npart <= 3L))

idx1 <- which(npart == 1L)
stopifnot(setequal(names(taeGut2)[idx1], seqlevels1))

idx2 <- which(npart == 2L)
m2 <- matrix(unlist(tmp[idx2]), ncol=2L, byrow=TRUE)
stopifnot(all(m2[ , 1L] == "chrUn"))
oo2 <- order(m2[ , 2L])
seqlevels2 <- names(taeGut2)[idx2[oo2]]

idx3 <- which(npart == 3L)
m3 <- matrix(unlist(tmp[idx3]), ncol=3L, byrow=TRUE)
stopifnot(all(m3[ , 3L] == "random"))
chrom <- factor(m3[ , 1L], levels=seqlevels1)
stopifnot(!anyNA(chrom))
oo3 <- order(chrom, m3[ , 2L])
seqlevels3 <- names(taeGut2)[idx3[oo3]]

taeGut2_seqlevels <- c(seqlevels1, seqlevels3, seqlevels2)

export(taeGut2[taeGut2_seqlevels], "taeGut2.sorted.2bit")

