make_test_GenomicFeature <- function() {
    new("GenomicFeature",
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
        strand = Rle(factor(strand(c("-", "+", "*", NA, "+", "-"))),
                     c(1, 2, 1, 1, 3, 2)),
        values = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

test_GenomicFeature_construction <- function() {
    checkException(GenomicFeature(letters), silent = TRUE)
    checkException(GenomicFeature(ranges = IRanges(1:10, 1:10)), silent = TRUE)
    checkException(GenomicFeature(letters, IRanges(1:10, 1:10)), silent = TRUE)
    checkException(GenomicFeature(letters, IRanges(1:26, 1:26),
                                 strand = letters), silent = TRUE)
    checkException(GenomicFeature(letters, IRanges(1:26, 1:26), score = 1:10),
                   silent = TRUE)
    checkException(GenomicFeature(letters, IRanges(1:26, 1:26), start = 1:26),
                   silent = TRUE)
    checkException(GenomicFeature(letters, IRanges(1:26, 1:26), end = 1:26),
                   silent = TRUE)
    checkException(GenomicFeature(letters, IRanges(1:26, 1:26), width = 1:26),
                   silent = TRUE)

    checkTrue(validObject(GenomicFeature()))
    checkTrue(validObject(GenomicFeature(letters, IRanges(1:26, 1:26))))
    checkTrue(validObject(GenomicFeature(letters, IRanges(1:26, 1:26),
                                        score = 1:26)))
    checkTrue(validObject(GenomicFeature(factor(letters), IRanges(1:26, 1:26))))
    checkTrue(validObject(GenomicFeature(1:10, IRanges(1:10, 1:10))))

    checkIdentical(GenomicFeature(seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
                                 ranges = IRanges(1:10, width = 10:1, names = head(letters,10)),
                                 strand = Rle(factor(strand(c("-", "+", "*", NA, "+", "-"))),
                                              c(1, 2, 1, 1, 3, 2)),
                                 score = 1:10, GC = seq(1, 0, length=10)),
                   make_test_GenomicFeature())
}

test_GenomicFeature_coercion <- function() {
    ## score, no strand
    gr <-
      GenomicFeature(seqnames = c(1,1,2),
                    ranges = IRanges(1:3,4:6, names = head(letters,3)),
                    score = c(10L,2L,NA))
    df <-
      data.frame(seqnames = as.character(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(rep(NA_character_, 3)),
                 score = c(10L,2L,NA),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(gr), df)

    ## strand, no score
    gr <-
      GenomicFeature(seqnames = c(1,1,2),
                    ranges = IRanges(1:3,4:6, names = head(letters,3)),
                    strand = strand(c("+", "-", "*")))
    df <-
      data.frame(seqnames = as.character(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(c("+", "-", "*")),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(gr), df)

    ## strand & score
    gr <-
      GenomicFeature(seqnames = c(1,1,2),
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
    checkIdentical(as.data.frame(gr), df)
}

test_GenomicFeature_accessors <- function() {
    ## seqnames
    checkException(seqnames(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GenomicFeature()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GenomicFeature()) <- letters,
                   silent = TRUE)

    gr <- make_test_GenomicFeature()
    val <- seqnames(gr)
    runValue(val) <- paste(runValue(val), ".new", sep="")
    seqnames(gr) <- val
    checkIdentical(seqnames(gr), val)

    gr <- make_test_GenomicFeature()
    val <- head(letters, length(gr))
    seqnames(gr) <- val
    checkIdentical(seqnames(gr), Rle(val))

    ## ranges
    checkException(ranges(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(ranges(make_test_GenomicFeature()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GenomicFeature()) <- IRanges(1:26, 1:26),
                   silent = TRUE)

    gr <- make_test_GenomicFeature()
    val <- IRanges(1:length(gr), width = 10)
    ranges(gr) <- val
    checkIdentical(ranges(gr), val)

    ## strand
    checkException(strand(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(strand(make_test_GenomicFeature()) <- NULL, silent = TRUE)
    checkException(strand(make_test_GenomicFeature()) <- letters, silent = TRUE)

    gr <- make_test_GenomicFeature()
    val <- Rle(strand("+"), length(gr))
    strand(gr) <- val
    checkIdentical(strand(gr), val)

    gr <- make_test_GenomicFeature()
    val <- rep(strand("+"), length(gr))
    strand(gr) <- val
    checkIdentical(strand(gr), Rle(val))

    ## values
    checkException(values(gr) <- DataFrame(strand = 1:length(gr)),
                   silent = TRUE)
    checkException(values(gr) <- DataFrame(score = letters), silent = TRUE)

    gr <- make_test_GenomicFeature()
    values(gr) <- NULL
    checkIdentical(values(gr),
                   new("DataFrame", nrows = length(gr), rownames = names(gr)))

    gr <- make_test_GenomicFeature()
    val <- DataFrame(x = 1:length(gr), y = head(letters, length(gr)))
    rownames(val) <- names(gr)
    values(gr) <- val
    checkTrue(validObject(gr))
    checkIdentical(values(gr), val)

    ## names
    checkException(names(gr) <- letters, silent = TRUE)

    gr <- make_test_GenomicFeature()
    names(gr) <- NULL
    checkIdentical(names(gr), NULL)

    gr <- make_test_GenomicFeature()
    names(gr) <- head(letters, length(gr))
    checkIdentical(names(gr), head(letters, length(gr)))
}

test_GenomicFeature_Ranges <- function() {
    ## start
    checkException(start(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(start(make_test_GenomicFeature()) <- letters, silent = TRUE)
    checkException(start(make_test_GenomicFeature()) <- 1:26, silent = TRUE)

    gr <- make_test_GenomicFeature()
    start(gr) <- as.numeric(seq_len(length(gr)))
    checkIdentical(start(gr), seq_len(length(gr)))

    ## end
    checkException(end(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(end(make_test_GenomicFeature()) <- letters, silent = TRUE)
    checkException(end(make_test_GenomicFeature()) <- 1:26, silent = TRUE)

    gr <- make_test_GenomicFeature()
    end(gr) <- as.numeric(10L + seq_len(length(gr)))
    checkIdentical(end(gr), 10L + seq_len(length(gr)))

    ## width
    checkException(width(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(width(make_test_GenomicFeature()) <- letters, silent = TRUE)
    checkException(width(make_test_GenomicFeature()) <- 1:26, silent = TRUE)

    gr <- make_test_GenomicFeature()
    width(gr) <- as.numeric(10L + seq_len(length(gr)))
    checkIdentical(width(gr), 10L + seq_len(length(gr)))
}

test_GenomicFeature_DataTable <- function() {
    checkIdentical(ncol(GenomicFeature()), 0L)
    checkIdentical(ncol(make_test_GenomicFeature()), 2L)

    checkException(colnames(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(colnames(make_test_GenomicFeature()) <- "a", silent = TRUE)
    checkException(colnames(make_test_GenomicFeature()) <- letters,
                   silent = TRUE)
    gr <- make_test_GenomicFeature()
    colnames(gr) <- c("a", "b")
    checkIdentical(colnames(gr), c("a", "b"))
}

test_GenomicFeature_Sequence <- function() {
    ## [
    gr <- make_test_GenomicFeature()
    checkException(gr[1000], silent = TRUE)
    checkException(gr["bad"], silent = TRUE)
    checkIdentical(gr, gr[])
    checkIdentical(as.data.frame(gr)[c(1,3,5),], as.data.frame(gr[c(1,3,5)]))
    checkIdentical(as.data.frame(gr)[c(1,3,5),-7],
                   as.data.frame(gr[c(1,3,5),"score"]))
    checkIdentical(as.data.frame(gr)[c(1,3,5),-7],
                   as.data.frame(gr[c(1,3,5),1]))

    ## [<-
    gr <- make_test_GenomicFeature()
    gr[] <- rev(gr)
    revgr <- rev(make_test_GenomicFeature())
    names(revgr) <- rev(names(revgr))
    checkIdentical(gr, revgr)

    ## c
    gr <- make_test_GenomicFeature()
    gr2 <- gr
    names(gr2) <- NULL
    checkException(c(GenomicFeature(), RangedData()), silent = TRUE)
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
    gr <- make_test_GenomicFeature()
    checkIdentical(gr[1:3], seqselect(gr, 1, 3))
    checkIdentical(gr[c(1:3, 1:3)], seqselect(gr, c(1,1), c(3,3)))

    ## seqselect<-
    gr1 <- make_test_GenomicFeature()
    gr1[1:3] <- make_test_GenomicFeature()[4:6]
    gr2 <- make_test_GenomicFeature()
    seqselect(gr2, 1, 3) <- make_test_GenomicFeature()[4:6]
    checkIdentical(gr1, gr2)

    ## window
    gr <- make_test_GenomicFeature()
    checkIdentical(gr[1:3], window(gr, 1, 3))
}
