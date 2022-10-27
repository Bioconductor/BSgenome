###
library(Biostrings)

### Download GCA_011100615.1_Macaca_fascicularis_6.0_genomic.fna.gz from
### https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/011/100/615/GCA_011100615.1_Macaca_fascicularis_6.0/
dna <- readDNAStringSet("GCA_011100615.1_Macaca_fascicularis_6.0_genomic.fna.gz")

### Check seqnames.
current_GenBankAccn <- unlist(heads(strsplit(names(dna), " ", fixed=TRUE), n=1L))
library(GenomeInfoDb)
chrominfo <- getChromInfoFromNCBI("GCA_011100615.1")
expected_GenBankAccn <- chrominfo[ , "GenBankAccn"]
stopifnot(setequal(expected_GenBankAccn, current_GenBankAccn))

### Reorder sequences.
dna <- dna[match(expected_GenBankAccn, current_GenBankAccn)]

### Rename sequences.
names(dna) <- chrominfo[ , "SequenceName"]

### Check sequence lengths.
expected_seqlengths <- chrominfo[ , "SequenceLength"]
stopifnot(all(width(dna) == expected_seqlengths | is.na(expected_seqlengths)))

### Export as 2bit.
library(rtracklayer)
export.2bit(dna, "Macaca_fascicularis_6.0.sorted.2bit")

