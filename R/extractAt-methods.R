### =========================================================================
### extractAt() methods
### -------------------------------------------------------------------------


.slow_extract_genome_sequences_at <- function(x, at)
{
    stopifnot(is(x, "BSgenome"), is(at, "GRanges"),
              identical(seqinfo(x), seqinfo(at)))

    seqlevels(at) <- seqlevelsInUse(at)
    grl <- split(at, seqnames(at))
    grl_seqlevels <- names(grl)

    ## Loop over the sequence names and extract the ranges.
    dnaset_list <- lapply(seq_along(grl),
        function(i) {
            gr <- grl[[i]]
            if (length(gr) == 0L)
                return(DNAStringSet())
            seqlevel <- grl_seqlevels[i]
            subject <- x[[seqlevel]]
            masks(subject) <- NULL
            is_circular <- isCircular(x)[[seqlevel]]
            ## 'seqlevel' will only be used for display in case of error.
            loadSubseqsFromStrandedSequence(subject, seqlevel,
                                            ranges(gr), strand(gr),
                                            is_circular)
        }
    )

    ## "unsplit" 'dnaset_list'.
    unsplit_list_of_XVectorList("DNAStringSet", dnaset_list,
                                as.factor(seqnames(at)))
}

### Use fast import() method for TwoBitFile objects.
.fast_extract_genome_sequences_at <- function(x, at)
{
    stopifnot(is(x, "BSgenome"), is(at, "GRanges"),
              identical(seqinfo(x), seqinfo(at)))

    ## import() doesn't seem safe. For example it accepts ranges with
    ## a start < 1L and returns a sequence with stuff in it that looks
    ## wrong! (but what would look right anyway?)
    ## So we check that the ranges in 'at' are valid, just to be safe.
    out_of_bound_idx <- GenomicRanges:::get_out_of_bound_index(at)
    if (length(out_of_bound_idx) != 0L)
        stop(wmsg("some ranges are out of bounds"))

    seqlevels_map <- setNames(names(x@user_seqnames), x@user_seqnames)
    seqlevels(at) <- seqlevels_map[seqlevels(at)]
    seqlevels(at) <- seqlevelsInUse(at)

    LRat <- splitLRgranges(at)
    ## import() ignores the strand of the ranges in 'which'.
    Lans <- import(x@single_sequences@twobitfile, which=LRat$L)
    Rans <- import(x@single_sequences@twobitfile, which=LRat$R)
    ans <- xscat(Lans, Rans)

    idx <- which(strand(at) == "-")
    ans[idx] <- reverseComplement(ans[idx])
    ans
}

### Must handle user defined seqlevels, injected SNPs, circularity, and strand.
.extract_genome_sequences_at <- function(x, at)
{
    stopifnot(is(x, "BSgenome"), is(at, "GRanges"))

    invalid_seqlevels <- setdiff(seqlevels(at), seqlevels(x))
    if (length(invalid_seqlevels) != 0L)
        stop(wmsg("'at' contains invalid chromosome name(s): ",
                  paste(invalid_seqlevels, sep=", ")))
    si <- merge(seqinfo(x), seqinfo(at))
    seqlevels(at) <- seqlevels(si)
    seqinfo(x) <- seqinfo(at) <- si

    if (is.null(SNPlocs_pkgname(x)) &&
        is(x@single_sequences, "TwobitNamedSequences"))
    {
        FUN <- .fast_extract_genome_sequences_at
    } else {
        FUN <- .slow_extract_genome_sequences_at
    }
    ans <- FUN(x, at)
    names(ans) <- names(at)
    ans
}

.extractAt_BSgenome <- function(x, at)
{
    if (is(at, "GRanges"))
        return(.extract_genome_sequences_at(x, at))
    if (is(at, "GRangesList")) {
        unlisted_at <- unlist(at, use.names=FALSE)
        unlisted_ans <- .extract_genome_sequences_at(x, unlisted_at)
        return(relist(unlisted_ans, at))
    }
    stop(wmsg("'at' must be a GRanges or GRangesList object"))
}

setMethod("extractAt", "BSgenome", .extractAt_BSgenome)

