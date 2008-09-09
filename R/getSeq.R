## If length(seqname) == 1:
##   o 'start' and 'end' are recycled to the length of the longest
##   o the result is a character vector (when 'as.XStringViews=FALSE')
##     or an XStringViews object (when 'as.XStringViews=TRUE')
##   o the length of the result is the length of the longest of 'start' and 'end'
##
## If length(seqname) != 1:
##   o 'start' and 'end' can't have more elements than 'seqname'
##   o 'start' and 'end' are recycled to the length of 'seqname' if necessary
##   o the result is a character vector ('as.XStringViews' arg is ignored)

getSeq <- function(bsgenome, seqname, start=NA, end=NA, as.XStringViews=FALSE)
{
    if (length(seqname) == 1) {
        ans <- Views(bsgenome[[seqname]], start=start, end=end)
        if (!as.XStringViews)
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
                                      as.XStringViews=FALSE))
    ans
}
