## Not exported
getDataAnnotationContribUrl <- function(type)
{
    ## Just to avoid codetools "no visible binding" NOTE
    biocinstallRepos <- NULL
    ## Sourcing this file will define the 'biocinstallRepos' function
    source("http://bioconductor.org/biocLite.R", local=TRUE)
    contrib.url(biocinstallRepos()[2], type=type)
}

available.genomes <- function(type=getOption("pkgType"))
{
    url <- getDataAnnotationContribUrl(type)
    pkgs <- available.packages(url)[, "Package"]
    pkgs <- pkgs[substr(pkgs, 1, 9) == "BSgenome."]
    names(pkgs) <- NULL
    pkgs
}

