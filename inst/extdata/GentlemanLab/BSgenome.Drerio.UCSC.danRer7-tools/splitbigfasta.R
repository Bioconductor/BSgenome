###
library(Biostrings)
danRer7 <- read.DNAStringSet("danRer7.fa")
idx1 <- grep("chr", names(danRer7), fixed=TRUE)
idx2 <- grep("Zv9_NA", names(danRer7), fixed=TRUE)
idx3 <- grep("Zv9_scaffold", names(danRer7), fixed=TRUE)

### Check that (idx1, idx2, idx3) forms a partition of seq_len(length(danRer7)).
stopifnot(identical(sort(c(idx1, idx2, idx3)), seq_len(length(danRer7))))

### Send each chromosome to a FASTA file.
seqnames <- paste("chr", c(1:25, "M"), sep="")
for (seqname in seqnames) {
    seq <- danRer7[match(seqname, names(danRer7))]
    filename <- paste(seqname, ".fa", sep="")
    cat("writing ", filename, "\n", sep="")
    write.XStringSet(seq, file=filename)
}

### Send all the Zv9_NA* sequences to a single FASTA file.
mseq2 <- danRer7[idx2]
mseq2 <- mseq2[order(as.integer(substr(names(mseq2), nchar("Zv9_NA")+1L, 999L)))]
write.XStringSet(mseq2, file="Zv9_NA.fa")

### Send all the Zv9_scaffold* sequences to a single FASTA file.
mseq3 <- danRer7[idx3]
mseq3 <- mseq3[order(as.integer(substr(names(mseq3), nchar("Zv9_scaffold")+1L, 999L)))]
write.XStringSet(mseq3, file="Zv9_scaffold.fa")

