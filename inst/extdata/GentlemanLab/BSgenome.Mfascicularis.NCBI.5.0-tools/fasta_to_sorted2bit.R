library(BSgenome)

INFILE <- "GCF_000364345.1_Macaca_fascicularis_5.0_genomic.fna.gz"
OUTFILE <- "Macaca_fascicularis_5.0.sorted.2bit"

assembly_report <- GenomeInfoDb:::fetch_assembly_report("GCF_000364345.1")
Mfascicularis <- readDNAStringSet(INFILE)

## Extract RefSeq accessions from long ugly names found in FASTA file.
refseq_accn <- as.character(phead(
                   CharacterList(
                       strsplit(names(Mfascicularis), " ", fixed=TRUE)
                   ),
                   n=1
               ))

## Replace long ugly names with official SequenceName.
m <- match(refseq_accn, assembly_report[ , "RefSeqAccn"])
stopifnot(all(!is.na(m)))
names(Mfascicularis) <- assembly_report[m, "SequenceName"]

## Sort Mfascicularis.
seqnames <- c(paste0("MFA", c(1:20, "X")), "MT")
scaffold_idx <- grep("^Scaffold", names(Mfascicularis))
scaffold_id <- sub("^Scaffold", "", names(Mfascicularis)[scaffold_idx])
seqnames <- c(seqnames, paste0("Scaffold", sort(as.integer(scaffold_id))))
Mfascicularis <- Mfascicularis[seqnames]

## Export as 2bit file.
export.2bit(Mfascicularis, OUTFILE)

