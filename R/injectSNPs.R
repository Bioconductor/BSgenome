### =========================================================================
### SNP injection
### -------------------------------------------------------------------------


available.SNPs <- function()
{
    url <- contrib.url(biocReposList()["aData"])
    pkgs <- available.packages(url)[, "Package"]
    pkgs <- pkgs[substr(pkgs, 1, 8) == "SNPlocs."]
    names(pkgs) <- NULL
    pkgs
}

setGeneric("injectSNPs", signature="x",
    function(x, SNPlocs_pkg) standardGeneric("injectSNPs")
)

setMethod("injectSNPs", "BSgenome",
    function(x, SNPlocs_pkg)
    {
        if (!is.null(SNPlocs_pkg(x)))
            stop("SNPs were already injected in genome 'x'. ",
                 "Injecting from more than 1 package is not supported.")
        if (!is.character(SNPlocs_pkg) || length(SNPlocs_pkg) != 1 || is.na(SNPlocs_pkg))
            stop("'SNPlocs_pkg' must be a single string")
        ans <- x
        ans@SNPlocs_pkg <- SNPlocs_pkg
        ans@.seqs_cache <- new.env(parent=emptyenv())
        ans@.link_counts <- new.env(parent=emptyenv())
        snp_count <- SNPcount(ans)
        if (!all(names(snp_count) %in% seqnames(ans)))
            stop("sequence names in package ", SNPlocs_pkg, " are not ",
                 "compatible with those in genome 'x'")
        ans
    }
)

