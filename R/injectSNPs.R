### =========================================================================
### SNP injection
### -------------------------------------------------------------------------

installed.SNPs <- function()
{
    pkgs <- installed.packages()[ , "Package"]
    pkgs <- pkgs[substr(pkgs, 1, 8) == "SNPlocs."]
    names(pkgs) <- NULL
    return(pkgs)
}

available.SNPs <- function(type=getOption("pkgType"))
{
    url <- getDataAnnotationContribUrl(type)
    pkgs <- available.packages(url)[, "Package"]
    pkgs <- pkgs[substr(pkgs, 1, 8) == "SNPlocs."]
    names(pkgs) <- NULL
    return(pkgs)
}

setGeneric("injectSNPs", signature="x",
    function(x, SNPlocs_pkgname) standardGeneric("injectSNPs")
)

setMethod("injectSNPs", "BSgenome",
    function(x, SNPlocs_pkgname)
    {
        if (!is.null(SNPlocs_pkgname(x)))
            stop("SNPs were already injected in genome 'x'. ",
                 "Injecting from more than 1 package is not supported.")
        if (!is.character(SNPlocs_pkgname) || length(SNPlocs_pkgname) != 1 || is.na(SNPlocs_pkgname))
            stop("'SNPlocs_pkgname' must be a single string")
        ans <- x
        ans@SNPlocs_pkgname <- SNPlocs_pkgname
        ans@.seqs_cache <- new.env(parent=emptyenv())
        ans@.link_counts <- new.env(parent=emptyenv())
        snp_count <- SNPcount(ans)
        if (!all(names(snp_count) %in% seqnames(ans)))
            stop("sequence names in package ", SNPlocs_pkgname, " are not ",
                 "compatible with those in genome 'x'")
        ans
    }
)

