### =========================================================================
### Some low-level internal (i.e. non-exported) utilities
### -------------------------------------------------------------------------


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### TODO: quickUnlist() and quickUnsplit() are of general interest and
### maybe should be moved somewhere else (to S4Vectors).
###
### Both functions assume that 'x' is a list of length >= 1 with no names,
### and that its list elements are of the same type.

quickUnlist <- function(x)
{
    x1 <- x[[1L]]
    if (is.factor(x1)) {
        ## Fast unlisting of a list of factors that have all
        ## the same levels in the same order.
        structure(unlist(x), class="factor", levels=levels(x1))
    } else {
        do.call(c, x)  # doesn't work on list of factors
    }
}

quickUnsplit <- function(x, f)
{
    idx <- split(seq_along(f), f)
    idx <- unlist(idx, use.names=FALSE)
    revidx <- integer(length(idx))
    revidx[idx] <- seq_along(idx)
    quickUnlist(x)[revidx]
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

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

