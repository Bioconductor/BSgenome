getSeq <- function(bsgenome, seqname, start=NA, end=NA)
{
    subBString(bsgenome[[seqname]], start, end)
}
