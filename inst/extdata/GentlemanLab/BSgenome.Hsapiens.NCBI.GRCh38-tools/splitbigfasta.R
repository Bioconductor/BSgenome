###
library(Biostrings)

### Download GCA_000001405.15_GRCh38_top-level.fna.gz from
### ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38
### then uncompress with gunzip GCA_000001405.15_GRCh38_top-level.fna.gz
GRCh38 <- readDNAStringSet("GCA_000001405.15_GRCh38_top-level.fna")
m_data <- strsplit(names(GRCh38), "|", fixed=TRUE)
fasta_headers <- matrix(unlist(m_data, use.names=FALSE), ncol=length(m_data))
accn <- fasta_headers[4L, ]
stopifnot(!any(duplicated(accn)))

assembly_report <- read.table("ftp://ftp.ncbi.nlm.nih.gov/genbank/genomes/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/GCA_000001405.15_GRCh38_assembly_report.txt", sep="\t", stringsAsFactors=FALSE)
accn2 <- assembly_report[[5L]]
stopifnot(!any(duplicated(accn2)))

stopifnot(setequal(accn, accn2))

seqnames <- assembly_report[[1L]][match(accn, accn2)]

for (i in seq_along(accn)) {
    filename <- paste0(seqnames[[i]], ".fa")
    cat("writing ", filename, "\n", sep="")
    writeXStringSet(GRCh38[i], file=filename, width=50L)
}

