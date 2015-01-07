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

