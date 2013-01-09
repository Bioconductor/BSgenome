###
library(Biostrings)
bosTau6 <- readDNAStringSet("bosTau6.fa")

### Partitioning:
is_chrUn <- grepl("^chrUn", names(bosTau6))
is_chrom <- !is_chrUn

### Send each chromosome to a FASTA file.
seqnames <- paste("chr", c(1:29, "X", "M"), sep="")
stopifnot(setequal(seqnames, names(bosTau6)[is_chrom]))
for (seqname in seqnames) {
    seq <- bosTau6[match(seqname, names(bosTau6))]
    filename <- paste(seqname, ".fa", sep="")
    cat("writing ", filename, "\n", sep="")
    writeXStringSet(seq, file=filename, width=50L)
}

### Send the 3286 chrUn_* sequences to 1 FASTA file.
chrUn_mseq <- bosTau6[is_chrUn]
writeXStringSet(chrUn_mseq, file="chrUn.fa", width=50L)

