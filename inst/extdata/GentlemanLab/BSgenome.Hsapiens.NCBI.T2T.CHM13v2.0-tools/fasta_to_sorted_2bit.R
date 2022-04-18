###
library(Biostrings)

### Download GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz from
### https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/009/914/755/GCA_009914755.4_T2T-CHM13v2.0/
dna <- readDNAStringSet("GCA_009914755.4_T2T-CHM13v2.0_genomic.fna.gz")

### Check seqnames.
current_GenBankAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L))
library(GenomeInfoDb)
chrominfo <- getChromInfoFromNCBI("GCA_009914755.4")
expected_GenBankAccn <- chrominfo[ , "GenBankAccn"]
stopifnot(setequal(expected_GenBankAccn, current_GenBankAccn))

### Reorder sequences.
dna <- dna[match(expected_GenBankAccn, current_GenBankAccn)]

### Rename sequences.
names(dna) <- chrominfo[ , "SequenceName"]

### Export as 2bit.
library(rtracklayer)
export.2bit(dna, "T2T-CHM13v2.0.sorted.2bit")

