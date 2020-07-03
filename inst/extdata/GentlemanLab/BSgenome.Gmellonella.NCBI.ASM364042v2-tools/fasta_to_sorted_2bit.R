###
library(Biostrings)

### Download GCF_003640425.2_ASM364042v2_genomic.fna.gz from
### https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/003/640/425/GCF_003640425.2_ASM364042v2/
dna <- readDNAStringSet("GCF_003640425.2_ASM364042v2_genomic.fna.gz")

### Check seqnames.
current_RefSeqAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L))
library(GenomeInfoDb)
chrominfo <- getChromInfoFromNCBI("GCF_003640425.2")
expected_RefSeqAccn <- chrominfo[ , "RefSeqAccn"]
stopifnot(setequal(expected_RefSeqAccn, current_RefSeqAccn))

### Reorder sequences.
dna <- dna[match(expected_RefSeqAccn, current_RefSeqAccn)]

### Rename sequences. An alternative would be to rename them to
### chrominfo[ , "SequenceName"] but these names are VERY ugly (e.g.
### "ScRZk8e_1;HRSCAF=1").
names(dna) <- expected_RefSeqAccn

### Export as 2bit.
library(rtracklayer)
export(dna, "ASM364042v2.sorted.2bit")

