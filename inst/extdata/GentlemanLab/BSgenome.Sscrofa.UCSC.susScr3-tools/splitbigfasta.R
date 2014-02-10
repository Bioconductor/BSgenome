###
library(Biostrings)
susScr3 <- readDNAStringSet("susScr3.fa")

### Partitioning:
is_chrom <- grepl("^chr", names(susScr3))
is_unplaced_scaffold <- grepl("^(GL|JH)", names(susScr3))

### Sanity check:
stopifnot(all(is_chrom | is_unplaced_scaffold))
stopifnot(!any(is_chrom & is_unplaced_scaffold))

### Send each chromosome to a FASTA file.
seqnames <- paste("chr", c(1:18, "X", "Y", "M"), sep="")
stopifnot(setequal(seqnames, names(susScr3)[is_chrom]))
for (seqname in seqnames) {
    seq <- susScr3[match(seqname, names(susScr3))]
    filename <- paste(seqname, ".fa", sep="")
    cat("writing ", filename, "\n", sep="")
    writeXStringSet(seq, file=filename, width=50L)
}

### Send the 4562 unplaced scaffolds to 1 FASTA file.
unplaced_scaffolds <- susScr3[is_unplaced_scaffold]
writeXStringSet(unplaced_scaffolds, file="unplaced_scaffolds.fa", width=50L)

