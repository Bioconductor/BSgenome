## If length(seqname) == 1:
##   o 'start' and 'end' are recycled to the length of the longest
##   o the result is a character vector (when 'as.BStringViews=FALSE')
##     or a BStringViews object (when 'as.BStringViews=TRUE')
##   o the length of the result is the length of the longest of 'start' and 'end'
##
## If length(seqname) != 1:
##   o 'start' and 'end' can't have more elements than 'seqname'
##   o 'start' and 'end' are recycled to the length of 'seqname' if necessary
##   o the result is a character vector ('as.BStringViews' arg is ignored)

getSeq <- function(bsgenome, seqname, start=NA, end=NA, as.BStringViews=FALSE)
{
    if (length(seqname) == 1) {
        ans <- views(bsgenome[[seqname]], start, end)
        if (!as.BStringViews)
            ans <- as.character(ans)
        return(ans)
    }
    lseqname <- length(seqname)
    # Adjust length(start)
    lstart <- length(start)
    if (lstart > lseqname)
        stop("'start' has more elements than 'seqname'")
    if (lstart < lseqname)
        start <- rep(start, length.out=lseqname)
    # Adjust length(end)
    lend <- length(end)
    if (lend > lseqname)
        stop("'end' has more elements than 'seqname'")
    if (lend < lseqname)
        end <- rep(end, length.out=lseqname)
    ans <- character(0)
    if (lseqname >= 1)
        for (i in 1:lseqname)
            ans <- append(ans, getSeq(bsgenome, seqname[i], start[i], end[i],
                                      as.BStringViews=FALSE))
    ans
}
