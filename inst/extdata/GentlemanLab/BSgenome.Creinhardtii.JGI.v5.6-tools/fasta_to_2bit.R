###
library(Biostrings)

### Download Creinhardtii_281_v5.0.fa.gz from
### https://phytozome-next.jgi.doe.gov/info/Creinhardtii_v5_6
dna <- readDNAStringSet("Creinhardtii_281_v5.0.fa.gz")

### Export as 2bit.
library(rtracklayer)
export.2bit(dna, "Creinhardtii_281_v5.0.2bit")

