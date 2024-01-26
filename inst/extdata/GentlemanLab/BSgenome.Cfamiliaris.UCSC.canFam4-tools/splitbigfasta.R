###
library(Biostrings)
canFam4 <- readDNAStringSet("canFam4.fa")

### Partitioning:
is_chrUn <- grepl("^chrUn", names(canFam4))
is_chrom <- !is_chrUn

### Send each chromosome to a FASTA file.
seqnames <- paste("chr", c(1:38, "X", "M"), sep="")
stopifnot(setequal(seqnames, names(canFam4)[is_chrom]))
for (seqname in seqnames) {
    seq <- canFam4[match(seqname, names(canFam4))]
    filename <- paste(seqname, ".fa", sep="")
    cat("writing ", filename, "\n", sep="")
    writeXStringSet(seq, file=filename, width=50L)
}

### Send the 1439 chrUn_* sequences to 1 FASTA file.
chrUn_mseq <- canFam4[is_chrUn]
writeXStringSet(chrUn_mseq, file="chrUn.fa", width=50L)

