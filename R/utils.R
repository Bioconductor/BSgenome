### =========================================================================
### Some low-level internal (i.e. non-exported) utilities
### -------------------------------------------------------------------------

getDataAnnotationContribUrl <- function(type=getOption("pkgType"))
{
    if (!suppressWarnings(suppressMessages(require(BiocInstaller,
                                                   quietly=TRUE)))) {
        ## Sourcing this file will install and load the BiocInstaller package.
        suppressWarnings(source("http://bioconductor.org/biocLite.R",
                                local=TRUE))
    }
    contrib.url(biocinstallRepos()["BioCann"], type=type)
}

