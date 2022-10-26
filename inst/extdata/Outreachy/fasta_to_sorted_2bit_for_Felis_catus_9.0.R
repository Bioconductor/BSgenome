###
library(Biostrings)

### Download GCF_000181335_Felis_catus_9.0_genomic.fna.gz from
### https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/181/335/GCF_000181335.3_Felis_catus_9.0/GCF_000181335.3_Felis_catus_9.0_genomic.fna.gz
dna <- readDNAStringSet("GCF_000181335.3_Felis_catus_9.0_genomic.fna.gz")
head(names(dna)) #inspect first 6 names on object

### Check seqnames.
current_RefSeqAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L)) 
library(GenomeInfoDb)
chrominfo <- getChromInfoFromNCBI("GCF_000181335.3")
expected_RefSeqAccn <- chrominfo[ , "RefSeqAccn"]
stopifnot(setequal(expected_RefSeqAccn, current_RefSeqAccn))

### Reorder sequences.
dna <- dna[match(expected_RefSeqAccn, current_RefSeqAccn)]

### Rename sequences.
names(dna) <- chrominfo[ , "SequenceName"]

### Export as 2bit.
library(rtracklayer)
export.2bit(dna, "Felis_catus.9.0.sorted.2bit")
