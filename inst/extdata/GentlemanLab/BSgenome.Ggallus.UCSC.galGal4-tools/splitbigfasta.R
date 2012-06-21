###
library(Biostrings)
galGal4 <- readDNAStringSet("galGal4.fa")

### Partitioning:
is_random <- grepl("^chr[^_]*_[^_]*_random$", names(galGal4))
is_chrUn <- grepl("^chrUn", names(galGal4))
is_chrom <- !(is_random | is_chrUn)

### Sanity check:
stopifnot(all(is_random | is_chrUn | is_chrom))
stopifnot(!any(is_random & is_chrUn))
stopifnot(!any(is_random & is_chrom))
stopifnot(!any(is_chrUn & is_chrom))

### Send each chromosome to a FASTA file.
seqnames <- paste("chr", c(1:28, 32, "M", "W", "Z", "LGE64", "LGE22C19W28_E50C23"), sep="")
stopifnot(setequal(seqnames, names(galGal4)[is_chrom]))
for (seqname in seqnames) {
    seq <- galGal4[match(seqname, names(galGal4))]
    filename <- paste(seqname, ".fa", sep="")
    cat("writing ", filename, "\n", sep="")
    writeXStringSet(seq, file=filename, width=50L)
}

### Send the 1805 chrNN_*_random sequences to 1 FASTA file.
random_mseq <- galGal4[is_random]
writeXStringSet(random_mseq, file="random.fa", width=50L)

### Send the 14093 chrUn_* sequences to 1 FASTA file.
chrUn_mseq <- galGal4[is_chrUn]
writeXStringSet(chrUn_mseq, file="chrUn.fa", width=50L)

