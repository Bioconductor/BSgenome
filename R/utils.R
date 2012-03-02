## Not exported
getDataAnnotationContribUrl <- function(type=getOption("pkgType"))
{
    ## Just to avoid codetools "no visible binding" NOTE
    biocinstallRepos <- NULL
    ## Sourcing this file will define the 'biocinstallRepos' function
    suppressWarnings(source("http://bioconductor.org/biocLite.R", local=TRUE))
    contrib.url(biocinstallRepos()["BioCann"], type=type)
}

