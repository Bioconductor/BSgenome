### =========================================================================
### available.genomes() and related utilities
### -------------------------------------------------------------------------
###


.BSGENOME_PREFIX <- "BSgenome."

.has_BSgenome_prefix <- function(x)
    substr(x, 1L, nchar(.BSGENOME_PREFIX)) == .BSGENOME_PREFIX


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### installed.genomes() and available.genomes()
###

.splitNameParts <- function(pkgs)
{
    if (length(pkgs) == 0L) {
        ans <- data.frame(pkgname=character(0),
                          organism=factor(),
                          provider=factor(),
                          provider_version=character(0),
                          masked=logical(0),
                          stringsAsFactors=FALSE)
        return(ans)
    }
    ## Remove ".masked" suffix.
    pkgs2 <- pkgs
    is_masked <- substr(pkgs2, nchar(pkgs2) - 6L, nchar(pkgs2)) == ".masked"
    idx <- which(is_masked)
    pkgs2[idx] <- substr(pkgs2[idx], 1L, nchar(pkgs2)[idx] - 7L)

    parts <- strsplit(pkgs2, ".", fixed=TRUE)
    nparts <- elementNROWS(parts)
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
    pkgs <- pkgs[.has_BSgenome_prefix(pkgs)]
    names(pkgs) <- NULL
    if (splitNameParts)
        pkgs <- .splitNameParts(pkgs)
    pkgs
}

available.genomes <- function(splitNameParts=FALSE, type=getOption("pkgType"))
{
    if (!isTRUEorFALSE(splitNameParts))
        stop("'splitNameParts' must be TRUE or FALSE")
    contrib_url <- get_data_annotation_contrib_url(type)
    pkgs <- available.packages(contrib_url)[, "Package"]
    pkgs <- pkgs[.has_BSgenome_prefix(pkgs)]
    names(pkgs) <- NULL
    if (splitNameParts)
        pkgs <- .splitNameParts(pkgs)
    pkgs
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### getBSgenome()
###

.stopOnAvailablePkg <- function(genome, is.source=FALSE)
    stop(genome, " package is not currently installed.\n",
         "  You first need to install it, which you can do with:\n",
         "      library(BiocInstaller)\n",
         "      biocLite(\"", genome, "\"",
         ifelse(is.source, ", type=\"source\"", ""), ")")

.stopOnMoreThanOneAvailablePkg <- function(genome, masked, is.source=FALSE)
    stop("Looks like there is more than one available ",
         ifelse(masked, "masked ", ""), "BSgenome data package\n",
         "  that matches genome assembly (a.k.a. provider version): ",
         genome, "\n",
         "  You first need to choose one (use 'available.genomes(",
         ifelse(is.source, "type=\"source\"", ""), ")' to get\n  the list). ",
         "Then install it, which you can do with:\n",
         "      library(BiocInstaller)\n",
         "      biocLite(\"<pkgname>\"",
         ifelse(is.source, ", type=\"source\"", ""), ")")

.getInstalledPkgnameFromProviderVersion <- function(genome, masked=FALSE)
{
    ## 1) Search installed packages.
    inst_pkgs <- installed.genomes(splitNameParts=TRUE)
    inst_pkgs <- inst_pkgs[inst_pkgs[ , "masked"] == masked, , drop=FALSE]
    idx <- which(genome == inst_pkgs[ , "provider_version"])
    if (length(idx) == 1L)
        return(inst_pkgs[idx , "pkgname"])
    if (length(idx) >= 2L)
        stop("Looks like you have more than one installed ",
             ifelse(masked, "masked ", ""), "BSgenome data package\n",
             "  that matches genome assembly (a.k.a. provider version): ",
             genome, "\n",
             "  Please disambiguate by specifying the full name of the ",
             "package you want\n  to use (use 'installed.genomes()' to ",
             "get the list).")

    ## 2) Search available packages.
    av_pkgs <- available.genomes(splitNameParts=TRUE)
    av_pkgs <- av_pkgs[av_pkgs[ , "masked"] == masked, ]
    idx <- which(genome == av_pkgs[ , "provider_version"])
    if (length(idx) == 1L) {
        genome <- av_pkgs[idx , "pkgname"]
        .stopOnAvailablePkg(genome)
    }
    if (length(idx) >= 2L)
        .stopOnMoreThanOneAvailablePkg(genome, masked)

    ## 3) Search available source packages.
    if (getOption("pkgType") != "source") {
        av_srcpkgs <- available.genomes(splitNameParts=TRUE, type="source")
        av_srcpkgs <- av_srcpkgs[av_srcpkgs[ , "masked"] == masked, ]
        idx <- which(genome == av_srcpkgs[ , "provider_version"])
        if (length(idx) == 1L) {
            genome <- av_srcpkgs[idx , "pkgname"]
            .stopOnAvailablePkg(genome, is.source=TRUE)
        }
        if (length(idx) >= 2L)
            .stopOnMoreThanOneAvailablePkg(genome, masked, is.source=TRUE)
    }

    ## All searches have failed.
    stop("Couldn't find a ",
         ifelse(masked, "masked ", ""), "BSgenome data package ",
         "that matches genome assembly\n  (a.k.a. provider version): ",
         genome, "\n\n",
         "  Please use 'available.genomes()' ",
         "(or 'available.genomes(type=\"source\")')\n",
         "  to check the list of BSgenome data packages that are available ",
         "in the\n  Bioconductor repositories for your version of ",
         "R/Bioconductor.\n  If you don't find what you are looking for, ",
         "please see the BSgenomeForge\n  vignette in the BSgenome software ",
         "package for how to forge a ",
         ifelse(masked, "masked ", ""), "BSgenome\n  data package ",
         "for your organism of interest.")
}

.getBSgenomeObjectFromInstalledPkgname <- function(genome)
{
    if (suppressWarnings(require(genome, quietly=TRUE, character.only=TRUE))) {
        ## 'genome' package is on the search list.
        pkgenvir <- as.environment(paste("package", genome, sep=":"))
        ans <- try(get(genome, envir=pkgenvir, inherits=FALSE), silent=TRUE)
        if (!is(ans, "BSgenome"))
            stop(genome, " doesn't look like a valid BSgenome data package")
        return(ans)
    }
    av_pkgs <- available.genomes()
    if (genome %in% av_pkgs)
        .stopOnAvailablePkg(genome)
    if (getOption("pkgType") != "source") {
        av_srcpkgs <- available.genomes(type="source")
        if (genome %in% av_srcpkgs)
            .stopOnAvailablePkg(genome, is.source=TRUE)
    }
    stop("Package ", genome, " is not available.\n",
         "  Please see the BSgenomeForge vignette in the BSgenome software ",
         "package\n  for how to forge a BSgenome data package for ",
         "your organism of interest.")
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
             "specifying a genome assembly (a.k.a.\n  provider version) ",
             "that refers unambiguously to an installed BSgenome data\n",
             "  package")
    if (!.has_BSgenome_prefix(genome)) {
        genome <- .getInstalledPkgnameFromProviderVersion(genome, masked=masked)
    } else if (masked) {
        warning("'masked' is ignored when 'genome' is supplied as ",
                "a full package name")
    }
    .getBSgenomeObjectFromInstalledPkgname(genome)
}

