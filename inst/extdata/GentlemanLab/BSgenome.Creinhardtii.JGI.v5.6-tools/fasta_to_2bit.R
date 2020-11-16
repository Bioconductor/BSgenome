###
library(Biostrings)

### Download Creinhardtii_281_v5.0.fa.gz from
### https://phytozome-next.jgi.doe.gov/info/Creinhardtii_v5_6
###   This data is public. Please cite the following:
###   Merchant, S. S., Prochnik, S. E., Vallon, O., Harris, E. H., Karpowicz, S. J., Witman, G. B., … Grossman, A. R. (2007). The Chlamydomonas Genome Reveals the Evolution of Key Animal and Plant Functions. Science, 318(5848), 245–250.
###   https://doi.org/10.1126/science.1143609
dna <- readDNAStringSet("Creinhardtii_281_v5.0.fa.gz")

### Export as 2bit.
library(rtracklayer)
export.2bit(dna, "Creinhardtii_281_v5.0.2bit")

