###
library(Biostrings)
rn5 <- readDNAStringSet("rn5.fa")

### Partitioning:
is_random <- grepl("^chr[^_]*_[^_]*_random$", names(rn5))
is_chrUn <- grepl("^chrUn", names(rn5))
is_chrom <- !(is_random | is_chrUn)

### Sanity check:
stopifnot(all(is_random | is_chrUn | is_chrom))
stopifnot(!any(is_random & is_chrUn))
stopifnot(!any(is_random & is_chrom))
stopifnot(!any(is_chrUn & is_chrom))

### Send each chromosome to a FASTA file.
seqnames <- paste("chr", c(1:20, "X", "M"), sep="")
stopifnot(setequal(seqnames, names(rn5)[is_chrom]))
for (seqname in seqnames) {
    seq <- rn5[match(seqname, names(rn5))]
    filename <- paste(seqname, ".fa", sep="")
    cat("writing ", filename, "\n", sep="")
    writeXStringSet(seq, file=filename, width=50L)
}

### Send the 1278 chrNN_*_random sequences to 1 FASTA file.
random_mseq <- rn5[is_random]
writeXStringSet(random_mseq, file="random.fa", width=50L)

### Send the 1439 chrUn_* sequences to 1 FASTA file.
chrUn_mseq <- rn5[is_chrUn]
writeXStringSet(chrUn_mseq, file="chrUn.fa", width=50L)

