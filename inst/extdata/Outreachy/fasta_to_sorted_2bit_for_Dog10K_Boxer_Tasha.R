###
library(Biostrings)

### Download GCF_000002285.5_Dog10K_Boxer_Tasha_genomic.fna.gz from
### https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002/285/GCF_000002285.5_Dog10K_Boxer_Tasha/
dna <- readDNAStringSet("GCF_000002285.5_Dog10K_Boxer_Tasha_genomic.fna.gz")

### Check seqnames.
current_RefSeqAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L))
library(GenomeInfoDb)
chrominfo <- getChromInfoFromNCBI("GCF_000002285.5")
expected_RefSeqAccn <- chrominfo[ , "RefSeqAccn"]
stopifnot(setequal(expected_RefSeqAccn, current_RefSeqAccn))

### Reorder sequences.
dna <- dna[match(expected_RefSeqAccn, current_RefSeqAccn)]

### Rename sequences.
names(dna) <- chrominfo[ , "SequenceName"]

### Export as 2bit.
library(rtracklayer)
export.2bit(dna, "Dog10K_Boxer_Tasha.sorted.2bit")
