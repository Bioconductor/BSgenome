.splitNameParts <- function(pkgs)
{
    tmp <- strsplit(pkgs, ".", fixed=TRUE)
    eltlens <- elementLengths(tmp)
    ## Some packages like BSgenome.Tgondii.ToxoDB.7.0 don't follow the rules
    ## and are made of more than 4 parts.
    has_more_than_4_parts <- which(eltlens > 4L)
    for (i in has_more_than_4_parts) {
        tmp_i <- tmp[[i]]
        tmp[[i]] <- c(tmp_i[1:3], paste(tmp_i[4:eltlens[i]], collapse="."))
    }
    ## And just in case one day we start to see some packages with less than
    ## 4 parts.
    has_less_than_4_parts <- which(eltlens < 4L)
    for (i in has_less_than_4_parts)
        tmp[[i]] <- c(tmp[[i]], rep.int(NA_character_, 4L - eltlens[i]))
    ## From now, any top-level element in 'tmp' is guaranteed to be a character
    ## vector of length 4.
    utmp <- unlist(tmp)
    idx1 <- (1:length(tmp)) * 4L - 3L
    cbind(pkgname=pkgs,
          organism=utmp[idx1+1L],
          provider=utmp[idx1+2L],
          provider_version=utmp[idx1+3L])
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

