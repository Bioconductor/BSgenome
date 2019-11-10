### =========================================================================
### Miscellaneous low-level utils
### -------------------------------------------------------------------------
###
### Unless stated otherwise, nothing in this file is exported.
###


get_data_annotation_contrib_url <- function(type=getOption("pkgType"))
{
    ## The BiocManager package is needed for repositories().
    if (!requireNamespace("BiocManager", quietly=TRUE)) {
        stop("Install 'BiocManager' from CRAN to get 'BioCann' contrib.url")
    }
    contrib.url(BiocManager::repositories()["BioCann"], type=type)
}

### TODO: Move this to GenomeInfoDb.
read_seqinfo_table <- function(filepath, genome=NA)
{
    df <- read.table(filepath, stringsAsFactors=FALSE)
    seqnames <- df[[1L]]
    Seqinfo(
        seqnames=seqnames,
        seqlengths=df[[2L]],
        isCircular=df[[3L]],
        genome=genome
    )
}

### Encode/decode character vector of single letters as/from raw object.
### Used in SNPlocs packages to store IUPAC ambiguity letters on disk as
### a raw vector with 1 byte per letter. Alternatively a DNAString object
### could be used for that!
encode_letters_as_bytes <- function(x) charToRaw(paste(x, collapse=""))

decode_bytes_as_letters <- function(x) rawToChar(x, multiple=TRUE)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Used in the context of sequence extraction from circular chromosomes
###

### Return 2 IRanges objects parallel to 'x'.
splitLRranges <- function(x, seqlength, seqname)
{
    stopifnot(is(x, "IntegerRanges"),
              isSingleInteger(seqlength),
              seqlength >= 1L)

    x_start <- start(x)
    x_end <- end(x)
    start0 <- x_start - 1L  # 0-based start
    shift <- start0 %% seqlength - start0
    end1 <- x_end + shift
    L_start <- x_start + shift
    L_end <- pmin(end1, seqlength)
    R_start <- rep.int(1L, length(x))
    R_end <- pmax(end1 - seqlength, 0L)
    idx <- which(R_end > seqlength)
    if (length(idx) != 0L)
        stop(wmsg("A region to extract from circular ",
                  "sequence \"", seqname, "\" crosses the end/start ",
                  "of the circle more than once. ",
                  "This is not supported at the moment, sorry!"))
    Lans <- IRanges(L_start, L_end)
    Rans <- IRanges(R_start, R_end)
    list(L=Lans, R=Rans)
}

### Return 2 GRanges objects parallel to 'x'.
splitLRgranges <- function(x)
{
    stopifnot(is(x, "GenomicRanges"))

    x_seqnames_id <- as.integer(seqnames(x))
    seqlengths <- unname(seqlengths(x))[x_seqnames_id]
    if (any(width(x) != 0L & seqlengths == 0L))
        stop(wmsg("cannot extract stuff from empty sequences"))
    x_is_circ <- unname(isCircular(x)) %in% TRUE
    split_idx <- which(x_is_circ[x_seqnames_id] & width(x) != 0L)

    Lans <- Rans <- granges(x)
    width(Rans) <- 0L
    if (length(split_idx) == 0L)
        return(list(L=Lans, R=Rans))

    x_start <- start(x)[split_idx]
    x_end <- end(x)[split_idx]
    seqlengths <- seqlengths[split_idx]
    start0 <- x_start - 1L  # 0-based start
    shift <- start0 %% seqlengths - start0
    end1 <- x_end + shift
    L_start <- x_start + shift
    L_end <- pmin(end1, seqlengths)
    R_start <- rep.int(1L, length(split_idx))
    R_end <- pmax(end1 - seqlengths, 0L)
    idx <- which(R_end > seqlengths)
    if (length(idx) != 0L)
        stop(wmsg("A region to extract from a circular ",
                  "sequence crosses the end/start ",
                  "of the circle more than once. ",
                  "This is not supported at the moment, sorry!"))
    ranges(Lans)[split_idx] <- IRanges(L_start, L_end)
    ranges(Rans)[split_idx] <- IRanges(R_start, R_end)
    list(L=Lans, R=Rans)
}

