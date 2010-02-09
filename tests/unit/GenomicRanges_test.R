make_test_GenomicRanges <- function() {
    new("GenomicRanges",
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
        strand = Rle(factor(strand(c("-", "+", "*", NA, "+", "-"))),
                     c(1, 2, 1, 1, 3, 2)),
        values = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

test_GenomicRanges_construction <- function() {
    checkException(GenomicRanges(letters), silent = TRUE)
    checkException(GenomicRanges(ranges = IRanges(1:10, 1:10)), silent = TRUE)
    checkException(GenomicRanges(letters, IRanges(1:10, 1:10)), silent = TRUE)
    checkException(GenomicRanges(letters, IRanges(1:26, 1:26),
                                 strand = letters), silent = TRUE)
    checkException(GenomicRanges(letters, IRanges(1:26, 1:26), score = 1:10),
                   silent = TRUE)
    checkException(GenomicRanges(letters, IRanges(1:26, 1:26), start = 1:26),
                   silent = TRUE)
    checkException(GenomicRanges(letters, IRanges(1:26, 1:26), end = 1:26),
                   silent = TRUE)
    checkException(GenomicRanges(letters, IRanges(1:26, 1:26), width = 1:26),
                   silent = TRUE)

    checkTrue(validObject(GenomicRanges()))
    checkTrue(validObject(GenomicRanges(letters, IRanges(1:26, 1:26))))
    checkTrue(validObject(GenomicRanges(letters, IRanges(1:26, 1:26),
                                        score = 1:26)))
    checkTrue(validObject(GenomicRanges(factor(letters), IRanges(1:26, 1:26))))
    checkTrue(validObject(GenomicRanges(1:10, IRanges(1:10, 1:10))))

    checkIdentical(GenomicRanges(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
                                 ranges = IRanges(1:10, width = 10:1, names = head(letters,10)),
                                 strand = Rle(factor(strand(c("-", "+", "*", NA, "+", "-"))),
                                              c(1, 2, 1, 1, 3, 2)),
                                 score = 1:10, GC = seq(1, 0, length=10)),
                   make_test_GenomicRanges())
}

test_GenomicRanges_coercion <- function() {
    ## score, no strand
    rd <- RangedData(IRanges(1:3,4:6), score = c(10L,2L,NA), space = c(1,1,2))
    rownames(rd) <- head(letters,3)
    gr <-
      GenomicRanges(seqnames = c(1,1,2),
                    ranges = IRanges(1:3,4:6, names = head(letters,3)),
                    score = c(10L,2L,NA))
    df <-
      data.frame(seqnames = as.character(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(rep(NA_character_, 3)),
                 score = c(10L,2L,NA),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as(rd, "GenomicRanges"), gr)
    checkIdentical(as.data.frame(gr), df)

    ## strand, no score
    rd <-
      RangedData(IRanges(1:3,4:6),
                 strand = strand(c("+", "-", "*")), 
                 space = c(1,1,2))
    rownames(rd) <- head(letters,3)
    gr <-
      GenomicRanges(seqnames = c(1,1,2),
                    ranges = IRanges(1:3,4:6, names = head(letters,3)),
                    strand = strand(c("+", "-", "*")))
    df <-
      data.frame(seqnames = as.character(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(c("+", "-", "*")),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as(rd, "GenomicRanges"), gr)
    checkIdentical(as.data.frame(gr), df)

    ## strand & score
    rd <-
      RangedData(IRanges(1:3,4:6),
                 strand = strand(c("+", "-", "*")), 
                 score = c(10L,2L,NA), 
                 space = c(1,1,2))
    rownames(rd) <- head(letters,3)
    gr <-
      GenomicRanges(seqnames = c(1,1,2),
                    ranges = IRanges(1:3,4:6, names = head(letters,3)),
                    strand = strand(c("+", "-", "*")),
                    score = c(10L,2L,NA))
    df <-
      data.frame(seqnames = as.character(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(c("+", "-", "*")),
                 score = c(10L,2L,NA),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as(rd, "GenomicRanges"), gr)
    checkIdentical(as.data.frame(gr), df)
}

test_GenomicRanges_accessors <- function() {
    ## seqnames
    checkException(seqnames(GenomicRanges()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GenomicRanges()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GenomicRanges()) <- letters,
                   silent = TRUE)

    gr <- make_test_GenomicRanges()
    val <- seqnames(gr)
    runValue(val) <- paste(runValue(val), ".new", sep="")
    seqnames(gr) <- val
    checkIdentical(seqnames(gr), val)

    gr <- make_test_GenomicRanges()
    val <- head(letters, length(gr))
    seqnames(gr) <- val
    checkIdentical(seqnames(gr), Rle(val))

    ## ranges
    checkException(ranges(GenomicRanges()) <- NULL, silent = TRUE)
    checkException(ranges(make_test_GenomicRanges()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GenomicRanges()) <- IRanges(1:26, 1:26),
                   silent = TRUE)

    gr <- make_test_GenomicRanges()
    val <- IRanges(1:length(gr), width = 10)
    ranges(gr) <- val
    checkIdentical(ranges(gr), val)

    ## strand
    checkException(strand(GenomicRanges()) <- NULL, silent = TRUE)
    checkException(strand(make_test_GenomicRanges()) <- NULL, silent = TRUE)
    checkException(strand(make_test_GenomicRanges()) <- letters, silent = TRUE)

    gr <- make_test_GenomicRanges()
    val <- Rle(strand("+"), length(gr))
    strand(gr) <- val
    checkIdentical(strand(gr), val)

    gr <- make_test_GenomicRanges()
    val <- rep(strand("+"), length(gr))
    strand(gr) <- val
    checkIdentical(strand(gr), Rle(val))

    ## values
    checkException(values(gr) <- DataFrame(strand = 1:length(gr)),
                   silent = TRUE)
    checkException(values(gr) <- DataFrame(score = letters), silent = TRUE)

    gr <- make_test_GenomicRanges()
    values(gr) <- NULL
    checkIdentical(values(gr),
                   new("DataFrame", nrows = length(gr), rownames = names(gr)))

    gr <- make_test_GenomicRanges()
    val <- DataFrame(x = 1:length(gr), y = head(letters, length(gr)))
    rownames(val) <- names(gr)
    values(gr) <- val
    checkTrue(validObject(gr))
    checkIdentical(values(gr), val)

    ## names
    checkException(names(gr) <- letters, silent = TRUE)

    gr <- make_test_GenomicRanges()
    names(gr) <- NULL
    checkIdentical(names(gr), NULL)

    gr <- make_test_GenomicRanges()
    names(gr) <- head(letters, length(gr))
    checkIdentical(names(gr), head(letters, length(gr)))
}

test_GenomicRanges_Ranges <- function() {
    ## start
    checkException(start(GenomicRanges()) <- NULL, silent = TRUE)
    checkException(start(make_test_GenomicRanges()) <- letters, silent = TRUE)
    checkException(start(make_test_GenomicRanges()) <- 1:26, silent = TRUE)

    gr <- make_test_GenomicRanges()
    start(gr) <- as.numeric(seq_len(length(gr)))
    checkIdentical(start(gr), seq_len(length(gr)))

    ## end
    checkException(end(GenomicRanges()) <- NULL, silent = TRUE)
    checkException(end(make_test_GenomicRanges()) <- letters, silent = TRUE)
    checkException(end(make_test_GenomicRanges()) <- 1:26, silent = TRUE)

    gr <- make_test_GenomicRanges()
    end(gr) <- as.numeric(10L + seq_len(length(gr)))
    checkIdentical(end(gr), 10L + seq_len(length(gr)))

    ## width
    checkException(width(GenomicRanges()) <- NULL, silent = TRUE)
    checkException(width(make_test_GenomicRanges()) <- letters, silent = TRUE)
    checkException(width(make_test_GenomicRanges()) <- 1:26, silent = TRUE)

    gr <- make_test_GenomicRanges()
    width(gr) <- as.numeric(10L + seq_len(length(gr)))
    checkIdentical(width(gr), 10L + seq_len(length(gr)))
}

test_GenomicRanges_DataTable <- function() {
    checkIdentical(ncol(GenomicRanges()), 0L)
    checkIdentical(ncol(make_test_GenomicRanges()), 2L)

    checkException(colnames(GenomicRanges()) <- NULL, silent = TRUE)
    checkException(colnames(make_test_GenomicRanges()) <- "a", silent = TRUE)
    checkException(colnames(make_test_GenomicRanges()) <- letters,
                   silent = TRUE)
    gr <- make_test_GenomicRanges()
    colnames(gr) <- c("a", "b")
    checkIdentical(colnames(gr), c("a", "b"))
}

test_GenomicRanges_Sequence <- function() {
    ## [
    gr <- make_test_GenomicRanges()
    checkException(gr[1000], silent = TRUE)
    checkException(gr["bad"], silent = TRUE)
    checkIdentical(gr, gr[])
    checkIdentical(as.data.frame(gr)[c(1,3,5),], as.data.frame(gr[c(1,3,5)]))
    checkIdentical(as.data.frame(gr)[c(1,3,5),-7],
                   as.data.frame(gr[c(1,3,5),"score"]))
    checkIdentical(as.data.frame(gr)[c(1,3,5),-7],
                   as.data.frame(gr[c(1,3,5),1]))

    ## [<-
    gr <- make_test_GenomicRanges()
    gr[] <- rev(gr)
    revgr <- rev(make_test_GenomicRanges())
    names(revgr) <- rev(names(revgr))
    checkIdentical(gr, revgr)

    ## c
    gr <- make_test_GenomicRanges()
    gr2 <- gr
    names(gr2) <- NULL
    checkException(c(GenomicRanges(), RangedData()), silent = TRUE)
    checkException(c(gr, gr[,-1]), silent = TRUE)
    checkIdentical(as.data.frame(c(gr, gr)),
                   as.data.frame(gr)[rep(seq_len(length(gr)), 2),])
    checkIdentical(as.data.frame(c(gr, gr2)),
                   rbind(as.data.frame(gr), as.data.frame(gr2)))
    checkIdentical(as.data.frame(c(gr2, gr)),
                   rbind(as.data.frame(gr2), as.data.frame(gr)))

    ## length
    checkIdentical(length(gr), length(gr@seqnames))

    ## seqselect
    gr <- make_test_GenomicRanges()
    checkIdentical(gr[1:3], seqselect(gr, 1, 3))
    checkIdentical(gr[c(1:3, 1:3)], seqselect(gr, c(1,1), c(3,3)))

    ## seqselect<-
    gr1 <- make_test_GenomicRanges()
    gr1[1:3] <- make_test_GenomicRanges()[4:6]
    gr2 <- make_test_GenomicRanges()
    seqselect(gr2, 1, 3) <- make_test_GenomicRanges()[4:6]
    checkIdentical(gr1, gr2)

    ## window
    gr <- make_test_GenomicRanges()
    checkIdentical(gr[1:3], window(gr, 1, 3))
}
