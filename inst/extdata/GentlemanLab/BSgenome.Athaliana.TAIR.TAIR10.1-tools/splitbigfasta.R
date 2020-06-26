###
library(Biostrings)
library(GenomeInfoDb)

### Download GCF_000001735.4_TAIR10.1_genomic.fna.gz from
### https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/735/GCF_000001735.4_TAIR10.1/
TAIR10.1 <- readDNAStringSet("GCF_000001735.4_TAIR10.1_genomic.fna.gz")
current_RefSeqAccn <- unlist(heads(strsplit(names(TAIR10.1), " ", fixed=TRUE), n=1L))
chrominfo <- getChromInfoFromNCBI("TAIR10.1")
expected_RefSeqAccn <- chrominfo[ , "RefSeqAccn"]
stopifnot(identical(expected_RefSeqAccn, current_RefSeqAccn))
names(TAIR10.1) <- chrominfo[ , "SequenceName"]

for (i in seq_along(TAIR10.1)) {
    filename <- paste0(names(TAIR10.1)[[i]], ".fa")
    cat("writing ", filename, "\n", sep="")
    writeXStringSet(TAIR10.1[i], file=filename, width=50L)
}

