make_test_GenomicFeatureList <- function() {
    GenomicFeatureList(
        a =
        new("GenomicFeature",
            seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
            ranges = IRanges(1:10, width = 10:1, names = head(letters, 10)),
            strand = Rle(strand(c("-", "+", "*", NA, "+", "-")), c(1, 2, 1, 1, 3, 2)),
            values = DataFrame(score = 1:10, GC = seq(1, 0, length=10))),
        b =
        new("GenomicFeature",
            seqnames = Rle(c("chr2", "chr4", "chr5"), c(3, 6, 4)),
            ranges = IRanges(1:13, width = 13:1, names = tail(letters, 13)),
            strand = Rle(strand(c("-", "+", "-")), c(4, 5, 4)),
            values = DataFrame(score = 1:13, GC = seq(0, 1, length=13))))
}

test_GenomicFeatureList_construction <- function() {
    checkException(GenomicFeatureList(IRangesList()), silent = TRUE)

    checkTrue(validObject(GenomicFeatureList()))
    checkTrue(validObject(GenomicFeatureList(GenomicFeature())))
    checkTrue(validObject(GenomicFeatureList(a = GenomicFeature())))
    checkTrue(validObject(make_test_GenomicFeatureList()))
}

test_GenomicFeatureList_coercion <- function() {
    ## RangedData -> GenomicFeatureList
    rd <-
      RangedData(space = c(1,1,2),
                 ranges = IRanges(1:3,4:6, names = head(letters,3)),
                 strand = strand(c("+", "-", "*")),
                 score = c(10L,2L,NA))
    gfl <-
      split(GenomicFeature(seqnames = c(1,1,2),
                           ranges = IRanges(1:3,4:6, names = head(letters,3)),
                           strand = strand(c("+", "-", "*")),
                           score = c(10L,2L,NA)), head(letters,3))
   checkIdentical(as(rd, "GenomicFeatureList"), gfl)

    ## as.data.frame
    gf1 <-
      GenomicFeature(seqnames = c(1,1,2),
                     ranges = IRanges(1:3,4:6, names = head(letters,3)),
                     strand = strand(c("+", "-", "*")),
                     score = c(10L,2L,NA))
    gf2 <-
      GenomicFeature(seqnames = c("chr1", "chr2"),
                     ranges = IRanges(1:2,1:2, names = tail(letters,2)),
                     score = 12:13)
    gfl <- GenomicFeatureList(a = gf1, b = gf2)
    df <-
      data.frame(feature = rep(c("a","b"), c(3, 2)),
                 seqnames = c(1,1,2,"chr1","chr2"),
                 start = c(1:3,1:2), end = c(4:6,1:2),
                 width = c(4L, 4L, 4L, 1L, 1L),
                 strand = strand(c("+", "-", "*", rep(NA_character_, 2))),
                 score = c(10L,2L,NA,12:13),
                 row.names = c(head(letters,3), tail(letters,2)),
                 stringsAsFactors = FALSE)
    checkIdentical(as.data.frame(gfl), df)
}

test_GenomicFeatureList_accessors <- function() {
    gfl <- make_test_GenomicFeatureList()
    checkIdentical(seqnames(gfl), RleList(lapply(gfl, seqnames), compress=TRUE))
    checkIdentical(ranges(gfl), IRangesList(lapply(gfl, ranges)))
    checkIdentical(strand(gfl), RleList(lapply(gfl, strand), compress=TRUE))
    checkIdentical(values(gfl), SplitDataFrameList(lapply(gfl, values)))
}

test_GenomicFeatureList_RangesList <- function() {
    gfl <- make_test_GenomicFeatureList()
    checkIdentical(start(gfl), IntegerList(lapply(gfl, start)))
    checkIdentical(end(gfl), IntegerList(lapply(gfl, end)))
    checkIdentical(width(gfl), IntegerList(lapply(gfl, width)))

    ## start
    checkException(start(GenomicFeatureList()) <- NULL, silent = TRUE)
    checkException(start(make_test_GenomicFeatureList()) <- 1:26, silent = TRUE)

    gfl <- make_test_GenomicFeatureList()
    orig <- start(gfl)
    start(gfl) <- orig + 1L
    checkIdentical(start(gfl), orig + 1L)

    ## end
    checkException(end(GenomicFeatureList()) <- NULL, silent = TRUE)
    checkException(end(make_test_GenomicFeatureList()) <- 1:26, silent = TRUE)

    gfl <- make_test_GenomicFeatureList()
    orig <- end(gfl)
    end(gfl) <- orig + 1L
    checkIdentical(end(gfl), orig + 1L)

    ## width
    checkException(width(GenomicFeatureList()) <- NULL, silent = TRUE)
    checkException(width(make_test_GenomicFeatureList()) <- 1:26, silent = TRUE)

    gfl <- make_test_GenomicFeatureList()
    orig <- width(gfl)
    width(gfl) <- orig + 1L
    checkIdentical(width(gfl), orig + 1L)
}

test_GenomicFeatureList_SplitDataFrameList <- function() {
    checkIdentical(ncol(GenomicFeatureList()), 0L)
    checkIdentical(ncol(make_test_GenomicFeatureList()), 2L)

    checkException(colnames(GenomicFeatureList()) <- NULL, silent = TRUE)
    checkException(colnames(make_test_GenomicFeatureList()) <- "a", silent = TRUE)
    checkException(colnames(make_test_GenomicFeatureList()) <- letters,
                   silent = TRUE)
    gfl <- make_test_GenomicFeatureList()
    colnames(gfl) <- c("a", "b")
    checkIdentical(colnames(gfl), c("a", "b"))

    gfl <- make_test_GenomicFeatureList()
    checkIdentical(gfl, gfl[])
    checkIdentical(gfl[,"score"],
                   GenomicFeatureList(lapply(gfl, function(x) x[,"score"])))
    checkIdentical(gfl[seqnames(gfl) == "chr2",],
                   GenomicFeatureList(lapply(gfl, function(x) 
                                      x[seqnames(x) == "chr2",])))
    checkIdentical(gfl[seqnames(gfl) == "chr2", "score"],
                   GenomicFeatureList(lapply(gfl, function(x) 
                                      x[seqnames(x) == "chr2", "score"])))
}
