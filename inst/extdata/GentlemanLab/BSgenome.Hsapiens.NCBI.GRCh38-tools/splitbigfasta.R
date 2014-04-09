###
library(Biostrings)
library(GenomeInfoDb)

### Download GCA_000001405.15_GRCh38_top-level.fna.gz from
### ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38
### then uncompress with gunzip GCA_000001405.15_GRCh38_top-level.fna.gz
GRCh38 <- readDNAStringSet("GCA_000001405.15_GRCh38_top-level.fna")
m_data <- strsplit(names(GRCh38), "|", fixed=TRUE)
fasta_headers <- matrix(unlist(m_data, use.names=FALSE), ncol=length(m_data))
GRCh38_accns <- fasta_headers[4L, ]
stopifnot(!any(duplicated(GRCh38_accns)))

GRCh38_accn2seqlevel <- GenomeInfoDb:::fetch_GenBankAccn2seqlevel_for_GRCh38()

stopifnot(setequal(GRCh38_accns, names(GRCh38_accn2seqlevel)))
GRCh38_seqlevels <- GRCh38_accn2seqlevel[GRCh38_accns]

for (i in seq_along(GRCh38)) {
    filename <- paste0(GRCh38_seqlevels[[i]], ".fa")
    cat("writing ", filename, "\n", sep="")
    writeXStringSet(GRCh38[i], file=filename, width=50L)
}

