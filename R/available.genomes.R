### =========================================================================
### available.genomes() and related utilities
### -------------------------------------------------------------------------
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### installed.genomes() and available.genomes()
###

.splitNameParts <- function(pkgs)
{
    ## Remove ".masked" suffix.
    pkgs2 <- pkgs
    is_masked <- substr(pkgs2, nchar(pkgs2) - 6L, nchar(pkgs2)) == ".masked"
    idx <- which(is_masked)
    pkgs2[idx] <- substr(pkgs2[idx], 1L, nchar(pkgs2)[idx] - 7L)

    parts <- strsplit(pkgs2, ".", fixed=TRUE)
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
               masked=is_masked,
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getBSgenome()
###

.stopIfPkgIsAvailable <- function(genome, is.source=FALSE)
    stop(genome, " package is not currently installed.\n",
         "  You first need to install it, which you can do with:\n",
         "      library(BiocInstaller)\n",
         "      biocLite(\"", genome, "\"",
         ifelse(is.source, ", type=\"source\"", ""), ")")

.stopIfMoreThanOnePkgIsAvailable <- function(genome, masked, is.source=FALSE)
    stop("Looks like there is more than one available ",
         ifelse(masked, "masked ", ""), "BSgenome data package\n",
         "  that matches genome assembly ", genome, ".\n",
         "  You first need to choose one (use 'available.genomes(",
         ifelse(is.source, "type=\"source\"", ""), ")' to get\n  the list). ",
         "Then install it, which you can do with:\n",
         "      library(BiocInstaller)\n",
         "      biocLite(\"<pkgname>\"",
         ifelse(is.source, ", type=\"source\"", ""), ")")

.stopWithHelpfulMsgAboutAvailablePkgs <- function(genome, masked=FALSE)
{
    av_pkgs <- available.genomes(splitNameParts=TRUE)
    ## 1) Try full match.
    if (genome %in% av_pkgs[ , "pkgname"])
        .stopIfPkgIsAvailable(genome)
    if (getOption("pkgType") != "source") {
        av_srcpkgs <- available.genomes(splitNameParts=TRUE, type="source")
        if (genome %in% av_srcpkgs[ , "pkgname"])
            .stopIfPkgIsAvailable(genome, is.source=TRUE)
    }
    ## 2) Try match on "provider_version".
    av_pkgs <- av_pkgs[av_pkgs[ , "masked"] == masked, ]
    idx <- which(genome == av_pkgs[ , "provider_version"])
    if (length(idx) == 1L) {
        genome <- av_pkgs[idx , "pkgname"]
        .stopIfPkgIsAvailable(genome)
    }
    if (length(idx) >= 2L)
        .stopIfMoreThanOnePkgIsAvailable(genome, masked)
    if (getOption("pkgType") != "source") {
        av_srcpkgs <- av_srcpkgs[av_srcpkgs[ , "masked"] == masked, ]
        idx <- which(genome == av_srcpkgs[ , "provider_version"])
        if (length(idx) == 1L) {
            genome <- av_srcpkgs[idx , "pkgname"]
            .stopIfPkgIsAvailable(genome, is.source=TRUE)
        }
        if (length(idx) >= 2L)
            .stopIfMoreThanOnePkgIsAvailable(genome, masked, is.source=TRUE)
    }
    ## 1) and 2) have failed.
    stop(ifelse(masked, "Masked ", ""), "BSgenome data package ", genome,
         " is not available.\n",
         "  See the BSgenomeForge vignette in the BSgenome software package ",
         "for how to\n  forge a BSgenome data package for this genome.")
}

.getBSgenomeObjectFromInstalledPkg <- function(genome, masked=FALSE)
{
    pkgenvir <- try(as.environment(paste("package", genome, sep=":")),
                    silent=TRUE)
    if (!is(pkgenvir, "try-error")) {
        ## 'genome' package is on the search list.
        if (masked)
            warning("'masked' is only used to disambiguate between valid ",
                    "BSgenome data packages\n  when 'genome' is specifying ",
                    "a genome assembly. Otherwise it's ignored.")
        ans <- try(get(genome, envir=pkgenvir, inherits=FALSE), silent=TRUE)
        if (!is(ans, "BSgenome"))
            stop(genome, " doesn't look like a valid BSgenome data package")
        return(ans)
    }

    ## Try to find an installed BSgenome data package that matches 'genome'.
    inst_pkgs <- installed.genomes(splitNameParts=TRUE)
    ##   1) Try full match.
    if (genome %in% inst_pkgs[ , "pkgname"]) {
        library(genome, character.only=TRUE)
        return(.getBSgenomeObjectFromInstalledPkg(genome, masked=masked))
    }
    ##   2) Try match on "provider_version".
    inst_pkgs <- inst_pkgs[inst_pkgs[ , "masked"] == masked, ]
    idx <- which(genome == inst_pkgs[ , "provider_version"])
    if (length(idx) == 1L) {
        genome <- inst_pkgs[idx , "pkgname"]
        return(.getBSgenomeObjectFromInstalledPkg(genome))
    }
    if (length(idx) >= 2L)
        stop("Looks like you have more than one installed ",
             ifelse(masked, "masked ", ""), "BSgenome data package\n",
             "  that matches genome assembly ", genome, ".\n",
             "  Please disambiguate by specifying the full name of the ",
             "package you want\n  to use (use 'installed.genomes()' to ",
             "get the list).")

    ## Try to find an available BSgenome data package that matches 'genome'
    ## and fail gracefully with a helpful error message.
    .stopWithHelpfulMsgAboutAvailablePkgs(genome, masked)
}

getBSgenome <- function(genome, masked=FALSE)
{
    if (!isTRUEorFALSE(masked))
        stop("'masked' must be TRUE or FALSE")
    if (is(genome, "BSgenome"))
        return(genome)
    if (!isSingleString(genome))
        stop("'genome' must be a BSgenome object, or the full name of an ",
             "installed\n  BSgenome data package, or a short string ",
             "specifying a genome assembly that\n  refers unambiguously ",
             "to an installed BSgenome data package")
    .getBSgenomeObjectFromInstalledPkg(genome, masked=masked)
}

