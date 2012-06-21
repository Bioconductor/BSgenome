###
library(Biostrings)
panTro3 <- read.DNAStringSet("panTro3.fa")

### Send each chromosome to a FASTA file.
seqnames <- paste("chr", c(1, "2A", "2B", 3:22, "X", "Y", "M"), sep="")
for (seqname in seqnames) {
    seq <- panTro3[match(seqname, names(panTro3))]
    filename <- paste(seqname, ".fa", sep="")
    cat("writing ", filename, "\n", sep="")
    write.XStringSet(seq, file=filename, width=50L)
}

### Send the chr*_random sequences to 1 FASTA file per sequence.
for (seqname in seqnames) {
    mseq <- panTro3[grep(paste0(seqname, "_"), names(panTro3), fixed=TRUE)]
    filename <- paste(seqname, "_random.fa", sep="")
    cat("writing ", filename, "\n", sep="")
    write.XStringSet(mseq, file=filename, width=50L)
}

### Send all the chrUn_* sequences to 1 FASTA file.
chrUn_mseq <- panTro3[grep("chrUn_", names(panTro3), fixed=TRUE)]
write.XStringSet(chrUn_mseq, file="chrUn.fa", width=50L)

