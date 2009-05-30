### =========================================================================
### getSeq()
### -------------------------------------------------------------------------
### TODO:
### - return extracted sequences in a DNAStringSet object;
### - some speed improvements (takes currently about 35 sec. to extract 1M
###   short sequences from BSgenome.Mmusculus.UCSC.mm9).

.getOneSeq <- function(bsgenome, name)
{
    if (name %in% seqnames(bsgenome))
        return(bsgenome[[name]])
    nhits <- 0L
    for (mseqname in mseqnames(bsgenome)) {
        mseq <- bsgenome[[mseqname]]
        ii <- grep(name, names(mseq))
        nhits <- nhits + length(ii)
        if (length(ii) == 1L)
            ans <- mseq[[ii]]
    }
    if (nhits == 0L)
        stop("sequence ", name, " not found")
    if (nhits > 1L)
        stop("sequence ", name, " found more than once, ",
             "please use a non-ambiguous name")
    ans
}

setGeneric("getSeq", function(x, ...) standardGeneric("getSeq"))

setMethod("getSeq", "BSgenome",
          function(x, names, start=NA, end=NA, width=NA, strand="+",
                   as.character=TRUE)
{
    if (missing(names)) {
        names <- seqnames(x)
    } else {
        if (is.factor(names))
            names <- as.vector(names)
        if (!is.character(names) || any(is.na(names)))
            stop("'names' must be a character vector (with no NAs)")
    }
    if (is.factor(strand))
        strand <- as.vector(strand)
    if (!is.character(strand) || !all(strand %in% c("+", "-")))
        stop("values in 'strand' must be \"+\"s or \"-\"s")
    if (!isTRUEorFALSE(as.character))
        stop("'as.character' must be TRUE or FALSE")
    l1 <- length(start)
    l2 <- length(end)
    l3 <- length(width)
    l4 <- length(strand)
    if (all(c(l1, l2, l3, l4) == 1L)) {
        if (length(names) == 0L) {
            if (as.character)
                ans <- character()
            else
                ans <- DNAStringSet()
            return(ans)
        }
        ans <- lapply(names, function(name)
                             subseq(.getOneSeq(x, name), start=start, end=end, width=width))
        ## length(ans) == length(names) >= 1L
        if (strand == "-")
            ans <- lapply(ans, reverseComplement)
        if (as.character)
            ## masks are removed before coercion to character vector
            return(sapply(ans, function(seq) {masks(seq) <- NULL; as.character(seq)}))
        if (length(ans) > 1L)
            stop("'as.character=FALSE' is not supported yet when extracting more than one sequence")
        ## length(ans) == length(names) == 1L
        return(ans[[1L]])
    }
    start <- IRanges:::.normargSEW(start, "start")
    end <- IRanges:::.normargSEW(end, "end")
    width <- IRanges:::.normargSEW(width, "width")
    l0 <- length(names)
    max01234 <- max(l0, l1, l2, l3, l4)
    if (max01234 == 0L) {
        if (as.character)
            ans <- character()
        else
            ans <- DNAStringSet()
        return(ans)
    }
    ## Recycling will fail for vectors of length 0
    if (l0 < max01234)
        names <- IRanges:::recycleVector(names, max01234)
    if (l1 < max01234)
        start <- IRanges:::recycleVector(start, max01234)
    if (l2 < max01234)
        end <- IRanges:::recycleVector(end, max01234)
    if (l3 < max01234)
        width <- IRanges:::recycleVector(width, max01234)
    if (l4 < max01234)
        strand <- IRanges:::recycleVector(strand, max01234)
    if (!as.character)
        stop("'as.character=FALSE' is not supported yet when extracting more than one sequence")
    ## The 4 lists belows have identical names (the REFSEQnames)
    REFSEQnames2start <- split(start, names, drop=TRUE)
    REFSEQnames2end <- split(end, names, drop=TRUE)
    REFSEQnames2width <- split(width, names, drop=TRUE)
    REFSEQnames2strand <- split(strand, names, drop=TRUE)
    REFSEQnames <- names(REFSEQnames2start)  # REFSEQnames has no duplicates
    if (!all(REFSEQnames %in% seqnames(x)))
        stop("the sequence names in 'names' must belong to 'seqnames(x)'")
    extractSeqsFromREFSEQ <- function(REFSEQname)
    {
        REFSEQ_length <- seqlengths(x)[REFSEQname]
        REFSEQ_start <- REFSEQnames2start[[REFSEQname]]
        REFSEQ_end <- REFSEQnames2end[[REFSEQname]]
        REFSEQ_width <- REFSEQnames2width[[REFSEQname]]
        REFSEQ_strand <- REFSEQnames2strand[[REFSEQname]]
        nseq <- length(REFSEQ_start)
        solved_SEW <- try(solveUserSEW(rep.int(REFSEQ_length, nseq),
                                       start=REFSEQ_start,
                                       end=REFSEQ_end,
                                       width=REFSEQ_width),
                          silent=TRUE)
        if (is(solved_SEW, "try-error"))
            stop("Invalid sequence coordinates.\n",
                 "  Please make sure the supplied 'start', 'end' and 'width' arguments\n",
                 "  are defining a region that is within the limits of the ", REFSEQname, " sequence.")
        subject <- x[[REFSEQname]]
        masks(subject) <- NULL
        ## All the views are guaranteed to be within the limits of the subject
        seqs <- as.character(Views(subject, solved_SEW))
        ## Very inefficient, sorry! (when we have an efficient "[<-" for
        ## XStringSet objects, we'll be able to do a much better job...)
        ii <- REFSEQ_strand == "-"
        seqs[ii] <- as.character(reverseComplement(DNAStringSet(seqs[ii])))
        seqs
    }
    REFSEQnames2seqs <- lapply(REFSEQnames, extractSeqsFromREFSEQ)
    unsplit(REFSEQnames2seqs, names, drop=TRUE)
})

