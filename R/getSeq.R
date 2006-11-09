getSeq <- function(bsgenome, seqname, start=NA, end=NA, as.character=TRUE)
{
    seq <- views(bsgenome[[seqname]], start, end)
    if (as.character)
        seq <- as.character(seq)
    seq
}
