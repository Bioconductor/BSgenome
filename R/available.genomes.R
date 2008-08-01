available.genomes <- function(type=getOption("pkgType"))
{
    url <- contrib.url(biocReposList()["aData"], type=type)
    pkgs <- available.packages(url)[, "Package"]
    pkgs <- pkgs[substr(pkgs, 1, 9) == "BSgenome."]
    names(pkgs) <- NULL
    pkgs
}

