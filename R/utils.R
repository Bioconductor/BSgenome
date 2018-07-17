### =========================================================================
### Some low-level internal (i.e. non-exported) utilities
### -------------------------------------------------------------------------


get_data_annotation_contrib_url <- function(type=getOption("pkgType"))
{
    ## The BiocManager package is needed for repositories().
    if (!requireNamespace("BiocManager", quietly=TRUE)) {
        ## Sourcing this file will install and load the BiocInstaller package.
        stop("Install 'BiocManager' from CRAN to get 'BioCann' contrib.url")
    }
    contrib.url(BiocManager::repositories()["BioCann"], type=type)
}

### TODO: Move this to GenomeInfoDb.
read_seqinfo_table <- function(filepath, genome=NA)
{
    df <- read.table(filepath, stringsAsFactors=FALSE)
    seqnames <- df[[1L]]
    Seqinfo(
        seqnames=seqnames,
        seqlengths=df[[2L]],
        isCircular=df[[3L]],
        genome=genome
    )
}

### Encode/decode character vector of single letters as/from raw object.
### Used in SNPlocs packages to store IUPAC ambiguity letters on disk as
### a raw vector with 1 byte per letter. Alternatively a DNAString object
### could be used for that!
encode_letters_as_bytes <- function(x) charToRaw(paste(x, collapse=""))

decode_bytes_as_letters <- function(x) rawToChar(x, multiple=TRUE)

