.splitNameParts <- function(pkgs)
{
    parts <- strsplit(pkgs, ".", fixed=TRUE)
    nparts <- elementLengths(parts)
    ## Some packages like BSgenome.Tgondii.ToxoDB.7.0 don't follow the rules
    ## and are made of more than 4 parts.
    has_more_than_4_parts <- which(nparts > 4L)
    for (i in has_more_than_4_parts) {
        tmp <- parts[[i]]
        parts[[i]] <- c(tmp[1:3], paste(tmp[4:length(tmp)], collapse="."))
    }
    ## And just in case one day we start to see some packages with less than
    ## 4 parts.
    has_less_than_4_parts <- which(nparts < 4L)
    for (i in has_less_than_4_parts) {
        tmp <- parts[[i]]
        parts[[i]] <- c(tmp, rep.int(NA_character_, 4L - length(tmp)))
    }
    ## From now, any top-level element in 'parts' is guaranteed to be a
    ## character vector of length 4.
    uparts <- unlist(parts)
    idx4 <- (1:length(parts)) * 4L
    data.frame(pkgname=pkgs,
               organism=factor(uparts[idx4 - 2L]),
               provider=factor(uparts[idx4 - 1L]),
               provider_version=uparts[idx4],
               stringsAsFactors=FALSE)
}

installed.genomes <- function(splitNameParts=FALSE)
{
    if (!isTRUEorFALSE(splitNameParts))
        stop("'splitNameParts' must be TRUE or FALSE")
    pkgs <- installed.packages()[ , "Package"]
    pkgs <- pkgs[substr(pkgs, 1, 9) == "BSgenome."]
    names(pkgs) <- NULL
    if (splitNameParts)
        pkgs <- .splitNameParts(pkgs)
    return(pkgs)
}

available.genomes <- function(splitNameParts=FALSE, type=getOption("pkgType"))
{
    if (!isTRUEorFALSE(splitNameParts))
        stop("'splitNameParts' must be TRUE or FALSE")
    url <- getDataAnnotationContribUrl(type)
    pkgs <- available.packages(url)[, "Package"]
    pkgs <- pkgs[substr(pkgs, 1, 9) == "BSgenome."]
    names(pkgs) <- NULL
    if (splitNameParts)
        pkgs <- .splitNameParts(pkgs)
    return(pkgs)
}

