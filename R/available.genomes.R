available.genomes <- function()
{
    url <- contrib.url(biocReposList()["aData"])
    pkgs <- available.packages(url)[, "Package"]
    pkgs <- pkgs[substr(pkgs, 1, 9) == "BSgenome."]
    names(pkgs) <- NULL
    pkgs
}

