###
library(Biostrings)
danRer6 <- read.DNAStringSet("danRer6.fa")
idx1 <- grep("chr", names(danRer6), fixed=TRUE)
idx2 <- grep("Zv8_NA", names(danRer6), fixed=TRUE)
idx3 <- grep("Zv8_scaffold", names(danRer6), fixed=TRUE)

### Check that (idx1, idx2, idx3) forms a partition of seq_len(length(danRer6)).
stopifnot(identical(sort(c(idx1, idx2, idx3)), seq_len(length(danRer6))))

### Send each chromosome to a FASTA file.
seqnames <- paste("chr", c(1:25, "M"), sep="")
for (seqname in seqnames) {
    seq <- danRer6[match(seqname, names(danRer6))]
    filename <- paste(seqname, ".fa", sep="")
    cat("writing ", filename, "\n", sep="")
    write.XStringSet(seq, file=filename)
}

### Send all the Zv8_NA* sequences to a single FASTA file.
mseq2 <- danRer6[idx2]
mseq2 <- mseq2[order(as.integer(substr(names(mseq2), nchar("Zv8_NA")+1L, 999L)))]
write.XStringSet(mseq2, file="Zv8_NA.fa")

### Send all the Zv8_scaffold* sequences to a single FASTA file.
mseq3 <- danRer6[idx3]
mseq3 <- mseq3[order(as.integer(substr(names(mseq3), nchar("Zv8_scaffold")+1L, 999L)))]
write.XStringSet(mseq3, file="Zv8_scaffold.fa")

