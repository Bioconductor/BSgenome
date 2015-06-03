library(BSgenome)

INFILE <- "GCF_000364345.1_Macaca_fascicularis_5.0_genomic.fna.gz"
OUTFILE <- "Macaca_fascicularis_5.0.sorted.2bit"

## Fetch assembly report from:
##   ftp://ftp.ncbi.nlm.nih.gov/genomes/ASSEMBLY_REPORTS/All/GCF_000364345.1.assembly.txt
assembly_report <- GenomeInfoDb:::fetch_assembly_report("GCF_000364345.1")
Mfascicularis <- readDNAStringSet(INFILE)

## Clean names on Mfascicularis to keep only the RefSeq accession.
names(Mfascicularis) <- as.character(phead(
                          CharacterList(
                            strsplit(names(Mfascicularis), " ", fixed=TRUE)
                          ),
                          n=1
                        ))

## Order sequences in Mfascicularis like in assembly report.
Mfascicularis <- Mfascicularis[assembly_report[ , "RefSeqAccn"]]

## Replace RefSeq accessions with official SequenceName from assembly report.
SequenceName <- assembly_report[ , "SequenceName"]
names(Mfascicularis) <- SequenceName

## Move MT sequence from last position to position after chromosomes (MFA*
## sequences) and before scaffolds (Scaffold* sequences).
oo <- c(grep("^MFA", SequenceName),
        grep("^MT$", SequenceName),
        grep("^Scaffold", SequenceName))
Mfascicularis <- Mfascicularis[oo]

## Export as 2bit file.
export.2bit(Mfascicularis, OUTFILE)

