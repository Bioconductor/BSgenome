###
library(Biostrings)

### Download GCF_000004195.4_UCB_Xtro_10.0_genomic.fna.gz from
### https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/004/195/GCF_000004195.4_UCB_Xtro_10.0/
dna <- readDNAStringSet("GCF_000004195.4_UCB_Xtro_10.0_genomic.fna.gz")

### Check seqnames.
current_RefSeqAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L))
library(GenomeInfoDb)
chrominfo <- getChromInfoFromNCBI("GCF_000004195.4")
expected_RefSeqAccn <- chrominfo[ , "RefSeqAccn"]
stopifnot(setequal(expected_RefSeqAccn, current_RefSeqAccn))

### Reorder sequences.
dna <- dna[match(expected_RefSeqAccn, current_RefSeqAccn)]

### Rename sequences.
names(dna) <- chrominfo[ , "SequenceName"]

### Export as 2bit.
library(rtracklayer)
export.2bit(dna, "UCB_Xtro_10.0.sorted.2bit")
