###
library(Biostrings)
rheMac3 <- readDNAStringSet("rheMac3.fa")

### Partitioning:
is_chrUn <- grepl("^chrUn", names(rheMac3))
is_chrom <- !is_chrUn

### Send each chromosome to a FASTA file.
seqnames <- paste("chr", c(1:20, "X", "M"), sep="")
stopifnot(setequal(seqnames, names(rheMac3)[is_chrom]))
for (seqname in seqnames) {
    seq <- rheMac3[match(seqname, names(rheMac3))]
    filename <- paste(seqname, ".fa", sep="")
    cat("writing ", filename, "\n", sep="")
    writeXStringSet(seq, file=filename, width=50L)
}

### Send the 34081 chrUn_* sequences to 1 FASTA file.
chrUn_mseq <- rheMac3[is_chrUn]
writeXStringSet(chrUn_mseq, file="chrUn.fa", width=50L)

