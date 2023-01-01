### =========================================================================
### SNPlocs objects
### -------------------------------------------------------------------------


setClass("SNPlocs",
    representation(
        "VIRTUAL",

        ## Provider of the SNPs (e.g. "dbSNP").
        provider="character",

        ## e.g. "dbSNP Human BUILD 149"
        provider_version="character",

        ## Official release date of the SNPs (e.g. "Nov 9, 2010").
        release_date="character",

        ## Official release name of the SNPs (e.g. "Build 132").
        release_name="character",

        ## URL to the place where the original SNP data was downloaded from.
        source_data_url="character",

        ## Date the original SNP data was downloaded.
        download_date="character",

        ## Reference genome of the SNPs.
        reference_genome="GenomeDescription",

        ## Named list of "sequence name translation tables" (one table per
        ## compatible genome and each table is represented by a named
        ## character vector).
        compatible_genomes="list",

        ## Name of the SNPlocs data package where the SNPlocs object is
        ## defined.
        data_pkgname="character"
    )
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Accessors
###

setMethod("provider", "SNPlocs", function(x) x@provider)

setMethod("providerVersion", "SNPlocs", function(x) x@provider_version)

setMethod("releaseDate", "SNPlocs", function(x) x@release_date)

setGeneric("releaseName", function(x) standardGeneric("releaseName"))

setMethod("releaseName", "SNPlocs", function(x) x@release_name)

setGeneric("referenceGenome", function(x) standardGeneric("referenceGenome"))

setMethod("referenceGenome", "SNPlocs", function(x) x@reference_genome)

setGeneric("compatibleGenomes",
    function(x) standardGeneric("compatibleGenomes")
)
setMethod("compatibleGenomes", "SNPlocs", function(x) x@compatible_genomes)

setMethod("organism", "SNPlocs",
    function(object) organism(referenceGenome(object))
)

setMethod("commonName", "SNPlocs",
    function(object) commonName(referenceGenome(object))
)

setMethod("seqinfo", "SNPlocs", function(x) seqinfo(referenceGenome(x)))

setMethod("seqnames", "SNPlocs", function(x) seqnames(referenceGenome(x)))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### The 'show' method
###

setMethod("show", "SNPlocs",
    function(object)
    {
        cat("# SNPlocs object for ", organism(object),
            " (", releaseName(object), ")\n", sep="")
        cat("# reference genome: ",
            providerVersion(referenceGenome(object)), "\n", sep="")
        cat("# nb of SNPs: ", sum(snpcount(object)), "\n", sep="")
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snpcount()
###

setGeneric("snpcount", function(x) standardGeneric("snpcount"))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### snplocs()
###
### Used internally for SNP injection. Not intended for the end user.
### Must return a 2-col data-frame-like object with columns "loc" (integer)
### and "alleles_as_ambig" (character).
###

setGeneric("snplocs", signature="x",
    function(x, seqname, ...) standardGeneric("snplocs")
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SNP extractors: snpsBySeqname(), snpsByOverlaps(), snpsById()
###

setGeneric("snpsBySeqname", signature="x",
    function(x, seqnames, ...) standardGeneric("snpsBySeqname")
)

### Same args and signature as GenomicFeatures::transcriptsByOverlaps()
### EXCEPT for 'minoverlap' default value that we set to zero so we also
### get SNPs that are insertions (relevant for XtraSNPlocs objects).
setGeneric("snpsByOverlaps", signature="x",
    function(x, ranges, ...) standardGeneric("snpsByOverlaps")
)

setGeneric("snpsById", signature="x",
    function(x, ids, ...) standardGeneric("snpsById")
)

### TODO: Avoid code duplication between normarg_ranges() and
### GenomicAlignments:::.normarg_param().
normarg_ranges <- function(ranges)
{
    if (isSingleString(ranges)) {
        tmp1 <- strsplit(ranges, ":", fixed=TRUE)[[1L]]
        if (length(tmp1) != 2L)
            stop(wmsg("when a character string, 'ranges' must be ",
                      "of the form \"ch14:5201-5300\""))
        tmp2 <- as.integer(strsplit(tmp1[2L], "-", fixed=TRUE)[[1L]])
        if (length(tmp2) != 2L || any(is.na(tmp2)))
            stop(wmsg("when a character string, 'ranges' must be ",
                      "of the form \"ch14:5201-5300\""))
        ranges <- GRanges(tmp1[1L], IRanges(tmp2[1L], tmp2[2L]))
        return(ranges)
    }
    if (!is(ranges, "GenomicRanges"))
        stop(wmsg("'ranges' ranges must be a GenomicRanges object ",
                  "or a character string of the form \"ch14:5201-5300\""))
    ranges
}

### Return an integer vector with no NAs parallel to 'ids'.
ids2rowids <- function(ids)
{
    if (!(is.character(ids) || is.numeric(ids)))
        stop(wmsg("'ids' must be a character or integer vector with no NAs"))
    if (S4Vectors:::anyMissing(ids))
        stop(wmsg("'ids' cannot contain NAs"))
    if (is.character(ids)) {
        prefixes <- unique(substr(ids, 1L, 2L))
        if ("rs" %in% prefixes) {
            if (!setequal(prefixes, "rs"))
                stop(wmsg("'ids' cannot mix SNP ids that are prefixed ",
                          "with \"rs\" with SNP ids that are not"))
            ## Drop the "rs" prefix.
            ids <- substr(ids, 3L, nchar(ids))
        }
        ids <- suppressWarnings(as.numeric(ids))
        if (S4Vectors:::anyMissing(ids))
            stop(wmsg("cannot extract the digital part of ",
                      "some SNP ids in 'ids'"))
    }
    if (length(ids) != 0L && min(ids) < 0)
        stop(wmsg("'ids' contains unrealistic SNP ids"))
    if (!is.integer(ids)) {
        ids <- suppressWarnings(as.integer(ids))
        if (S4Vectors:::anyMissing(ids))
            stop(wmsg("'ids' contains SNP ids that are too big"))
    }
    ids
}

### Return a list of 2 vectors parallel to each other. The 1st and 2nd vectors
### are respectively the row indices (integer vector) and the 'user_ids' vector
### (which can be character, numeric, or integer), with the not found entries
### removed from both of them.
### Note that, if 'ifnotfound="error"' then the 2 returned vectors are parallel
### to input vectors 'user_rowids' and 'user_ids'.
rowids2rowidx <- function(user_rowids, user_ids, x_rowids_env, ifnotfound)
{
    if (length(ls(x_rowids_env)) == 0L)
        stop(wmsg("BSgenome internal error: data contains no SNP ids"))
    rowidx <- lookup_rowids(user_rowids, x_rowids_env)
    notfound_idx <- which(is.na(rowidx))
    if (length(notfound_idx) != 0L) {
        if (length(notfound_idx) <= 10L) {
            ids_to_show <- user_ids[notfound_idx]
        } else {
            ids_to_show <- c(user_ids[notfound_idx[1:9]], "...")
        }
        ids_to_show <- paste0(ids_to_show, collapse=", ")
        if (ifnotfound == "error")
            stop(wmsg("SNP ids not found: ",
                      ids_to_show,
                      "\n\nUse 'ifnotfound=\"drop\"' to drop them."))
        if (ifnotfound == "warning")
            warning(wmsg("SNP ids not found: ", ids_to_show,
                         "\n\nThey were dropped."))
        rowidx <- rowidx[-notfound_idx]
        user_ids <- user_ids[-notfound_idx]
    }
    list(rowidx, user_ids)
}

