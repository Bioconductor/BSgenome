available.SNPs <- function()
{
    url <- contrib.url(biocReposList()["aData"])
    pkgs <- available.packages(url)[, "Package"]
    pkgs <- pkgs[substr(pkgs, 1, 8) == "SNPlocs."]
    names(pkgs) <- NULL
    pkgs
}

