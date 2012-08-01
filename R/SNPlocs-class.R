### =========================================================================
### SNPlocs objects
### -------------------------------------------------------------------------

setClass("SNPlocs",
    representation(
        ## provider of the SNPs (e.g. "dbSNP")
        provider="character",

        ## creation date (in compact format) of the flat files found at
        ## download_url (look at 2nd line of each file e.g.
        ## CREATED ON: 2012-06-08 14:53, and use the most recent date in
        ## case of mixed dates, e.g. "20120608")
        provider_version="character",

        ## official release date of the SNPs (e.g. "Nov 9, 2010")
        release_date="character",

        ## official release name of the SNPs (e.g. "Build 132")
        release_name="character",

        ## URL to the place where the original SNP data was downloaded
        download_url="character",

        ## date the original SNP data was downloaded
        download_date="character",

        ## reference genome of the SNPs
        reference_genome="GenomeDescription",

        ## named list of "sequence name translation tables" (one table per
        ## compatible genome and each table is represented by a named character
        ## vector)
        compatible_genomes="list",

        ## package name and absolute path to local directory where to find
        ## the serialized objects containing the SNPs
        data_pkgname="character",
        data_dirpath="character",
        data_serialized_objnames="character",

        .data_cache="environment"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

### Not intended to be used directly.
newSNPlocs <- function(provider, provider_version,
                       release_date, release_name,
                       download_url, download_date,
                       reference_genome, compatible_genomes,
                       data_pkgname, data_dirpath)
{
    data_serialized_objnames <- c(
        "SNPcount",
        "all_rsids",
        paste(seqlevels(reference_genome), "_snplocs", sep="")
    )
    new("SNPlocs",
        provider=provider,
        provider_version=provider_version,
        release_date=release_date,
        release_name=release_name,
        download_url=download_url,
        download_date=download_date,
        reference_genome=reference_genome,
        compatible_genomes=compatible_genomes,
        data_pkgname=data_pkgname,
        data_dirpath=data_dirpath,
        data_serialized_objnames=data_serialized_objnames,
        .data_cache=new.env(hash=TRUE, parent=emptyenv()))
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("provider", "SNPlocs", function(x) x@provider)

setMethod("providerVersion", "SNPlocs", function(x) x@provider_version)

setMethod("releaseDate", "SNPlocs", function(x) x@release_date)

setMethod("releaseName", "SNPlocs", function(x) x@release_name)

setGeneric("referenceGenome", function(x) standardGeneric("referenceGenome"))
setMethod("referenceGenome", "SNPlocs", function(x) x@reference_genome)

setGeneric("compatibleGenomes",
    function(x) standardGeneric("compatibleGenomes")
)
setMethod("compatibleGenomes", "SNPlocs", function(x) x@compatible_genomes)

setMethod("organism", "SNPlocs", function(x) organism(referenceGenome(x)))

setMethod("species", "SNPlocs", function(x) species(referenceGenome(x)))

setMethod("seqinfo", "SNPlocs", function(x) seqinfo(referenceGenome(x)))

setMethod("seqnames", "SNPlocs", function(x) seqnames(referenceGenome(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

setMethod("show", "SNPlocs",
    function(object)
    {
        cat(class(object), " object for ", organism(object), " (",
            provider(object), " ", releaseName(object), ")\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpcount() and snplocs()
###

### Creates objects "on-the-fly" (not serialized).
### WARNING: Improper calls to .get_SNPlocs_data() by the .create_object()
### function can lead to infinite recursive loops!
.create_object <- function(objname)
{
    if (objname == "empty_snplocs") {
        obj <- data.frame(loc=integer(0), alleles=raw(0))
        return(obj)
    }
    if (objname == "empty_ufsnplocs") {
        obj <- data.frame(RefSNP_id=character(0),
                          alleles_as_ambig=character(0),
                          loc=integer(0),
                          stringsAsFactors=FALSE)
        return(obj)
    }

    # add more here...
    stop("don't know how to create object '", objname, "'")
}

.get_SNPlocs_data <- function(x, objname, caching=TRUE)
{
    datacache <- x@.data_cache
    not_cached <- !exists(objname, envir=datacache)
    if (not_cached) {
        if (objname %in% x@data_serialized_objnames) {
            filename <- paste(objname, ".rda", sep="")
            filepath <- file.path(x@data_dirpath, filename)
            load(filepath, envir=datacache)
        } else {
            assign(objname, .create_object(objname), envir=datacache)
        }
    }
    ans <- get(objname, envir=datacache)
    if (not_cached && !caching)
        rm(list=objname, envir=datacache)
    ans
}

setGeneric("snpcount", function(x) standardGeneric("snpcount"))
setMethod("snpcount", "SNPlocs",
    function(x)
    {
        objname <- "SNPcount"
        ans <- .get_SNPlocs_data(x, objname)
        if (!is.integer(ans) || !identical(names(ans), seqlevels(x)))
            stop("internal error: '", objname, "' data set is broken.\n",
                 "       Please contact the maintainer of the ",
                 x@data_pkgname, "\n       package.")
        ans        
    }
)

### Load raw snplocs.
.load_raw_snplocs <- function(x, seqname, caching)
{
    objname <- paste(seqname, "_snplocs", sep="")
    ans <- .get_SNPlocs_data(x, objname, caching=caching)
    empty_snplocs <- .get_SNPlocs_data(x, "empty_snplocs")
    if (!identical(sapply(ans, class), sapply(empty_snplocs, class)))
        stop("internal error: unexpected col names and/or col types\n",
             "       for the '", objname, "' data set.\n",
             "       Please contact the maintainer of the ",
             x@data_pkgname, "\n       package.")
    if (nrow(ans) != snpcount(x)[seqname])
        stop("internal error: nb of rows in object '", objname, " ",
             "doesn't match the\n       nb of SNPs reported by snpcount().\n",
             "       Please contact the maintainer of the ",
             x@data_pkgname, "\n       package.")
    ans
}

.get_rsid_offsets <- function(x)
{
    offsets <- c(0L, cumsum(snpcount(x)))
    offsets <- offsets[-length(offsets)]
    names(offsets) <- names(snpcount(x))
    offsets
}

### Load rs ids for a given sequence. Return them in an integer vector.
.load_rsids <- function(x, seqname)
{
    all_rsids <- .get_SNPlocs_data(x, "all_rsids")
    seq_pos <- match(seqname, names(snpcount(x)))
    offset <- .get_rsid_offsets(x)[seq_pos]
    idx <- seq_len(snpcount(x)[seq_pos]) + offset
    all_rsids[idx]
}

### Get user-friendly snplocs.
.get_ufsnplocs <- function(x, seqname, caching)
{
    rsids <- as.character(.load_rsids(x, seqname))
    snplocs <- .load_raw_snplocs(x, seqname, caching)
    alleles <- safeExplode(rawToChar(snplocs$alleles))
    data.frame(RefSNP_id=rsids,
               alleles_as_ambig=alleles,
               loc=snplocs$loc,
               stringsAsFactors=FALSE)
}

.SNPlocsAsGranges <- function(x, ufsnplocs, seqname)
{
    if (is(seqname, "Rle")) {
        if (length(seqname) != nrow(ufsnplocs)
         || !identical(levels(seqname), seqlevels(x)))
            stop("when an Rle, 'seqname' must be a factor Rle ",
                 "of length 'nrow(ufsnplocs)' and levels 'SEQNAMES'")
        ans_seqnames <- seqname
    } else {
        if (!is.factor(seqname))
            seqname <- factor(seqname, levels=seqlevels(x))
        if (length(seqname) == 1L)
            ans_seqnames <- Rle(seqname, nrow(ufsnplocs))
        else if (length(seqname) == nrow(ufsnplocs))
            ans_seqnames <- Rle(seqname)
        else
            stop("'length(seqname)' must be 1 or 'nrow(ufsnplocs)'")
    }
    if (nrow(ufsnplocs) == 0L)
        ans_ranges <- IRanges()
    else
        ans_ranges <- IRanges(start=ufsnplocs$loc, width=1L)
    ans_strand <- Rle(strand("+"), nrow(ufsnplocs))
    ans <- GRanges(seqnames=ans_seqnames,
                   ranges=ans_ranges,
                   strand=ans_strand,
                   RefSNP_id=ufsnplocs$RefSNP_id,
                   alleles_as_ambig=ufsnplocs$alleles_as_ambig)
    seqinfo(ans) <- seqinfo(x)
    ans
}

setGeneric("snplocs", signature="x",
    function(x, seqname, ...) standardGeneric("snplocs")
)

### Returns a data frame (when 'as.GRanges=FALSE') or a GRanges object
### (when 'as.GRanges=TRUE').
setMethod("snplocs", "SNPlocs",
    function(x, seqname, as.GRanges=FALSE, caching=TRUE)
    {
        if (!is.character(seqname) || any(is.na(seqname)))
            stop("'seqname' must be a character vector with no NAs")
        if (!all(seqname %in% seqlevels(x)))
            stop("all 'seqname' elements must be in: ",
                 paste(seqlevels(x), collapse=", "))
        if (!isTRUEorFALSE(as.GRanges))
            stop("'as.GRanges' must be TRUE or FALSE")
        if (!isTRUEorFALSE(caching))
            stop("'caching' must be TRUE or FALSE")
        if (!as.GRanges) {
            if (length(seqname) != 1L)
                stop("'seqname' must be of length 1 when 'as.GRanges' is FALSE")
            return(.get_ufsnplocs(x, seqname, caching))
        }
        if (length(seqname) == 0L) {
            empty_ufsnplocs <- .get_SNPlocs_data(x, "empty_ufsnplocs")
            ans <- .SNPlocsAsGranges(x, empty_ufsnplocs, character(0))
            return(ans)
        }
        list_of_ufsnplocs <- lapply(seqname,
                                    .get_ufsnplocs, x=x, caching=caching)
        ufsnplocs <- do.call(rbind, list_of_ufsnplocs)
        seqnames <- Rle(factor(seqname, levels=seqlevels(x)),
                        unlist(lapply(list_of_ufsnplocs, nrow), use.names=FALSE))
        .SNPlocsAsGranges(x, ufsnplocs, seqnames)
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpid2loc(), snpid2alleles(), and snpid2grange()
###

.normargSnpid <- function(snpid)
{
    if (!is.vector(snpid))
        stop("'snpid' must be an integer or character vector")
    if (IRanges:::anyMissing(snpid))
        stop("'snpid' cannot contain NAs")
    if (is.numeric(snpid)) {
        if (!is.integer(snpid))
            snpid <- as.integer(snpid)
        return(snpid)
    }
    if (!is.character(snpid))
        stop("'snpid' must be an integer or character vector")
    prefixes <- unique(substr(snpid, 1L, 2L))
    if ("rs" %in% prefixes) {
        if (!setequal(prefixes, "rs"))
            stop("'snpid' cannot mix ids that are prefixed with \"rs\" ",
                 "with ids that are not")
        ## Drop the "rs" prefix.
        snpid <- substr(snpid, 3, nchar(snpid))
    }
    snpid <- suppressWarnings(as.integer(snpid))
    if (IRanges:::anyMissing(snpid))
        stop("cannot extract the digital part of some ids in 'snpid'")
    snpid
}

### Returns a named integer vector where each (name, value) pair corresponds
### to a supplied SNP id (typically an rs id). The name is the chromosome of
### the SNP id and the value is the row index in the serialized snplocs data
### frame corresponding to the SNP id.
.snpid2rowidx <- function(x, snpid)
{
    if (length(snpid) == 0L) {
        idx <- integer(0)
    } else {
        all_rsids <- .get_SNPlocs_data(x, "all_rsids")
        idx <- match(snpid, all_rsids)
        bad_snpid_idx <- which(is.na(idx))
        if (length(bad_snpid_idx) != 0L) {
            bad_snpid <- snpid[bad_snpid_idx]
            bad_snpid_in1string <- paste(bad_snpid, collapse=", ")
            stop("SNP id(s) not found: ", bad_snpid_in1string)
        }
    }
    seqidx <- findInterval(idx - 1L, cumsum(snpcount(x))) + 1L
    rowidx <- idx - .get_rsid_offsets(x)[seqidx]
    names(rowidx) <- names(snpcount(x))[seqidx]
    rowidx
}

### Returns a named integer vector where each (name, value) pair corresponds
### to a supplied SNP id (typically an rs id). The name is the chromosome of
### the SNP id and the value is its position on the chromosome.
setGeneric("snpid2loc", signature="x",
    function(x, snpid, ...) standardGeneric("snpid2loc")
)

setMethod("snpid2loc", "SNPlocs",
    function(x, snpid, caching=TRUE)
    {
        snpid <- .normargSnpid(snpid)
        if (!isTRUEorFALSE(caching))
            stop("'caching' must be TRUE or FALSE")
        rowidx <- .snpid2rowidx(x, snpid)
        if (length(rowidx) == 0L) {
            ans <- integer(0)
        } else {
            rowidx_list <- split(unname(rowidx), names(rowidx))
            loc_list <- lapply(names(rowidx_list),
                function(seqname) {
                    idx <- rowidx_list[[seqname]]
                    .load_raw_snplocs(x, seqname, caching)$loc[idx]
                })
            ans <- unsplit(loc_list, names(rowidx))
        }
        names(ans) <- names(rowidx)
        ans
    }
)

### Returns a named character vector where each (name, value) pair corresponds
### to a supplied SNP id (typically an rs id). The name is the chromosome of
### the SNP id and the value is a single IUPAC code representing the associated
### alleles.
setGeneric("snpid2alleles", signature="x",
    function(x, snpid, ...) standardGeneric("snpid2alleles")
)

setMethod("snpid2alleles", "SNPlocs",
    function(x, snpid, caching=TRUE)
    {
        snpid <- .normargSnpid(snpid)
        if (!isTRUEorFALSE(caching))
            stop("'caching' must be TRUE or FALSE")
        rowidx <- .snpid2rowidx(x, snpid)
        if (length(rowidx) == 0L) {
            ans <- raw(0)
        } else {
            rowidx_list <- split(unname(rowidx), names(rowidx))
            alleles_list <- lapply(names(rowidx_list),
                function(seqname) {
                    idx <- rowidx_list[[seqname]]
                    .load_raw_snplocs(x, seqname, caching)$alleles[idx]
                })
            ans <- unsplit(alleles_list, names(rowidx))
        }
        ans <- safeExplode(rawToChar(ans))
        names(ans) <- names(rowidx)
        ans
    }
)

setGeneric("snpid2grange", signature="x",
    function(x, snpid, ...) standardGeneric("snpid2grange")
)

setMethod("snpid2grange", "SNPlocs",
    function(x, snpid, caching=TRUE)
    {
        snpid <- .normargSnpid(snpid)
        if (!isTRUEorFALSE(caching))
            stop("'caching' must be TRUE or FALSE")
        loc <- snpid2loc(x, snpid, caching=caching)
        alleles <- snpid2alleles(x, snpid, caching=caching)
        ufsnplocs <- data.frame(RefSNP_id=as.character(snpid),
                                alleles_as_ambig=unname(alleles),
                                loc=unname(loc),
                                stringsAsFactors=FALSE)
        .SNPlocsAsGranges(x, ufsnplocs, names(loc))
    }
)

