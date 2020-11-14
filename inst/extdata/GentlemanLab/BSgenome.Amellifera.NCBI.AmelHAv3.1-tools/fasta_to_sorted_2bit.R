###
library(Biostrings)

### Download GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz from
### https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/254/395/GCF_003254395.2_Amel_HAv3.1/
dna <- readDNAStringSet("GCF_003254395.2_Amel_HAv3.1_genomic.fna.gz")

### Check seqnames.
current_RefSeqAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L))
library(GenomeInfoDb)
chrominfo <- getChromInfoFromNCBI("GCF_003254395.2")
expected_RefSeqAccn <- chrominfo[ , "RefSeqAccn"]
stopifnot(setequal(expected_RefSeqAccn, current_RefSeqAccn))

### Reorder sequences.
dna <- dna[match(expected_RefSeqAccn, current_RefSeqAccn)]

### Rename sequences.
names(dna) <- chrominfo[ , "SequenceName"]

### Export as 2bit.
library(rtracklayer)
export.2bit(dna, "Amel_HAv3.1.sorted.2bit")

