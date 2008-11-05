.getOneSeq <- function(bsgenome, name)
{
    if (name %in% seqnames(bsgenome))
        return(bsgenome[[name]])
    nhits <- 0L
    for (mseqname in mseqnames(bsgenome)) {
        mseq <- bsgenome[[mseqname]]
        ii <- grep(name, names(mseq))
        nhits <- nhits + length(ii)
        if (length(ii) == 1)
            ans <- mseq[[ii]]
    }
    if (nhits == 0)
        stop("sequence ", name, " not found")
    if (nhits > 1)
        stop("sequence ", name, " found more than once, ",
             "please use a non-ambiguous name")
    ans
}

getSeq <- function(bsgenome, names, start=NA, end=NA, width=NA, as.character=TRUE)
{
    if (!is(bsgenome, "BSgenome"))
        stop("'bsgenome' must be a BSgenome object")
    if (missing(names))
        names <- seqnames(bsgenome)
    else if (!is.character(names) || any(is.na(names)))
        stop("'names' must be a character vector (with no NAs)")
    if (length(names) == 0) {
        ans <- character(0)
        if (!as.character)
            ans <- DNAStringSet(ans)
        return(ans)
    }
    ans <- lapply(names, function(name)
                    subseq(.getOneSeq(bsgenome, name), start=start, end=end, width=width))
    ## length(ans) == length(names) >= 1
    if (as.character)
        ## masks are removed before coercion to character vector
        return(sapply(ans, function(seq) {masks(seq) <- NULL; as.character(seq)}))
    if (length(names) > 1)
        stop("'as.character=FALSE' is not supported when 'length(names) > 1'")
    ## length(ans) == length(names) == 1
    ans[[1]]
}

