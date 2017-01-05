### =========================================================================
### Some low-level internal (i.e. non-exported) utilities
### -------------------------------------------------------------------------


get_data_annotation_contrib_url <- function(type=getOption("pkgType"))
{
    ## The BiocInstaller package is needed for biocinstallRepos().
    if (!requireNamespace("BiocInstaller", quietly=TRUE)) {
        ## Sourcing this file will install and load the BiocInstaller package.
        suppressWarnings(source("http://bioconductor.org/biocLite.R",
                                local=TRUE))
    }
    contrib.url(BiocInstaller::biocinstallRepos()["BioCann"], type=type)
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

