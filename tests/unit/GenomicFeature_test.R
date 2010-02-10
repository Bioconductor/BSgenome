make_test_GenomicFeature <- function() {
    new("GenomicFeature",
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
        strand = Rle(strand(c("-", "+", "*", NA, "+", "-")), c(1, 2, 1, 1, 3, 2)),
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
    checkException(GenomicFeature(letters, IRanges(1:26, 1:26),
                                  feature = letters), silent = TRUE)

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
    ## no strand or score
    gf <-
      GenomicFeature(seqnames = c(1,1,2),
                     ranges = IRanges(1:3,4:6, names = head(letters,3)))
    df <-
      data.frame(seqnames = as.character(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(rep(NA_character_, 3)),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(gf), df)

    ## score, no strand
    gf <-
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
    checkIdentical(as.data.frame(gf), df)

    ## strand, no score
    gf <-
      GenomicFeature(seqnames = c(1,1,2),
                    ranges = IRanges(1:3,4:6, names = head(letters,3)),
                    strand = strand(c("+", "-", "*")))
    df <-
      data.frame(seqnames = as.character(c(1,1,2)),
                 start = 1:3, end = 4:6, width = c(4L, 4L, 4L),
                 strand = strand(c("+", "-", "*")),
                 row.names = head(letters,3),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(gf), df)

    ## strand & score
    gf <-
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
    checkIdentical(as.data.frame(gf), df)
}

test_GenomicFeature_accessors <- function() {
    ## seqnames
    checkException(seqnames(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GenomicFeature()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GenomicFeature()) <- letters,
                   silent = TRUE)

    gf <- make_test_GenomicFeature()
    val <- seqnames(gf)
    runValue(val) <- paste(runValue(val), ".new", sep="")
    seqnames(gf) <- val
    checkIdentical(seqnames(gf), val)

    gf <- make_test_GenomicFeature()
    val <- head(letters, length(gf))
    seqnames(gf) <- val
    checkIdentical(seqnames(gf), Rle(val))

    ## ranges
    checkException(ranges(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(ranges(make_test_GenomicFeature()) <- NULL, silent = TRUE)
    checkException(seqnames(make_test_GenomicFeature()) <- IRanges(1:26, 1:26),
                   silent = TRUE)

    gf <- make_test_GenomicFeature()
    val <- IRanges(1:length(gf), width = 10)
    ranges(gf) <- val
    checkIdentical(ranges(gf), val)

    ## strand
    checkException(strand(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(strand(make_test_GenomicFeature()) <- NULL, silent = TRUE)
    checkException(strand(make_test_GenomicFeature()) <- letters, silent = TRUE)

    gf <- make_test_GenomicFeature()
    val <- Rle(strand("+"), length(gf))
    strand(gf) <- val
    checkIdentical(strand(gf), val)

    gf <- make_test_GenomicFeature()
    val <- rep(strand("+"), length(gf))
    strand(gf) <- val
    checkIdentical(strand(gf), Rle(val))

    ## values
    checkException(values(gf) <- DataFrame(strand = 1:length(gf)),
                   silent = TRUE)
    checkException(values(gf) <- DataFrame(score = letters), silent = TRUE)

    gf <- make_test_GenomicFeature()
    values(gf) <- NULL
    checkIdentical(values(gf),
                   new("DataFrame", nrows = length(gf), rownames = names(gf)))

    gf <- make_test_GenomicFeature()
    val <- DataFrame(x = 1:length(gf), y = head(letters, length(gf)))
    rownames(val) <- names(gf)
    values(gf) <- val
    checkTrue(validObject(gf))
    checkIdentical(values(gf), val)

    ## names
    checkException(names(gf) <- letters, silent = TRUE)

    gf <- make_test_GenomicFeature()
    names(gf) <- NULL
    checkIdentical(names(gf), NULL)

    gf <- make_test_GenomicFeature()
    names(gf) <- head(letters, length(gf))
    checkIdentical(names(gf), head(letters, length(gf)))
}

test_GenomicFeature_Ranges <- function() {
    ## start
    checkException(start(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(start(make_test_GenomicFeature()) <- letters, silent = TRUE)
    checkException(start(make_test_GenomicFeature()) <- 1:26, silent = TRUE)

    gf <- make_test_GenomicFeature()
    start(gf) <- as.numeric(seq_len(length(gf)))
    checkIdentical(start(gf), seq_len(length(gf)))

    ## end
    checkException(end(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(end(make_test_GenomicFeature()) <- letters, silent = TRUE)
    checkException(end(make_test_GenomicFeature()) <- 1:26, silent = TRUE)

    gf <- make_test_GenomicFeature()
    end(gf) <- as.numeric(10L + seq_len(length(gf)))
    checkIdentical(end(gf), 10L + seq_len(length(gf)))

    ## width
    checkException(width(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(width(make_test_GenomicFeature()) <- letters, silent = TRUE)
    checkException(width(make_test_GenomicFeature()) <- 1:26, silent = TRUE)

    gf <- make_test_GenomicFeature()
    width(gf) <- as.numeric(10L + seq_len(length(gf)))
    checkIdentical(width(gf), 10L + seq_len(length(gf)))
}

test_GenomicFeature_DataTable <- function() {
    checkIdentical(ncol(GenomicFeature()), 0L)
    checkIdentical(ncol(make_test_GenomicFeature()), 2L)

    checkException(colnames(GenomicFeature()) <- NULL, silent = TRUE)
    checkException(colnames(make_test_GenomicFeature()) <- "a", silent = TRUE)
    checkException(colnames(make_test_GenomicFeature()) <- letters,
                   silent = TRUE)
    gf <- make_test_GenomicFeature()
    colnames(gf) <- c("a", "b")
    checkIdentical(colnames(gf), c("a", "b"))
}

test_GenomicFeature_Sequence <- function() {
    ## [
    gf <- make_test_GenomicFeature()
    checkException(gf[1000], silent = TRUE)
    checkException(gf["bad"], silent = TRUE)
    checkIdentical(gf, gf[])
    checkIdentical(as.data.frame(gf)[c(1,3,5),], as.data.frame(gf[c(1,3,5)]))
    checkIdentical(as.data.frame(gf)[c(1,3,5),-7],
                   as.data.frame(gf[c(1,3,5),"score"]))
    checkIdentical(as.data.frame(gf)[c(1,3,5),-7],
                   as.data.frame(gf[c(1,3,5),1]))

    ## [<-
    gf <- make_test_GenomicFeature()
    gf[] <- rev(gf)
    revgf <- rev(make_test_GenomicFeature())
    names(revgf) <- rev(names(revgf))
    checkIdentical(gf, revgf)

    ## c
    gf <- make_test_GenomicFeature()
    gf2 <- gf
    names(gf2) <- NULL
    checkException(c(GenomicFeature(), RangedData()), silent = TRUE)
    checkException(c(gf, gf[,-1]), silent = TRUE)
    checkIdentical(as.data.frame(c(gf, gf)),
                   as.data.frame(gf)[rep(seq_len(length(gf)), 2),])
    checkIdentical(as.data.frame(c(gf, gf2)),
                   rbind(as.data.frame(gf), as.data.frame(gf2)))
    checkIdentical(as.data.frame(c(gf2, gf)),
                   rbind(as.data.frame(gf2), as.data.frame(gf)))

    ## length
    checkIdentical(length(gf), length(gf@seqnames))

    ## seqselect
    gf <- make_test_GenomicFeature()
    checkIdentical(gf[1:3], seqselect(gf, 1, 3))
    checkIdentical(gf[c(1:3, 1:3)], seqselect(gf, c(1,1), c(3,3)))

    ## seqselect<-
    gf1 <- make_test_GenomicFeature()
    gf1[1:3] <- make_test_GenomicFeature()[4:6]
    gf2 <- make_test_GenomicFeature()
    seqselect(gf2, 1, 3) <- make_test_GenomicFeature()[4:6]
    checkIdentical(gf1, gf2)

    ## split
    gf <- make_test_GenomicFeature()
    checkException(split(gf, NULL), silent = TRUE)
    checkIdentical(split(unname(gf)),
                   GenomicFeatureList(lapply(structure(seq_len(length(gf)),
                                names = as.character(seq_len(length(gf)))),
                                             function(i) unname(gf)[i])))
    checkIdentical(split(gf),
                   GenomicFeatureList(lapply(structure(seq_len(length(gf)),
                                                       names = names(gf)),
                                             function(i) gf[i])))
    checkIdentical(split(gf, rep(c("a", "b"), each=5)),
                   GenomicFeatureList(a = head(gf, 5), b = tail(gf, 5)))

    ## window
    gf <- make_test_GenomicFeature()
    checkIdentical(gf[1:3], window(gf, 1, 3))
}
