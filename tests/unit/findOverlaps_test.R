test_matching_seqnames <- function()
{
    prefixes <- c("", "chr", "Chr", "CHR", "chrom", "Chrom", "CHROM")
    P <- length(prefixes)
    for (isArabic1 in c(TRUE, FALSE)) {
        for (pre1 in seq_len(P - 1L)) {
            for (pre2 in ((pre1 + 1L):P)) {
                num1 <- if(isArabic1) 1:24 else as.character(as.roman(1:30))
                num2 <- if(!isArabic1) as.character(as.roman(1:24)) else 1:30
                seq1 <- paste(prefixes[pre1], num1, sep = "")
                seq2 <- paste(prefixes[pre2], num2, sep = "")
                checkTrue(!BSgenome:::.similarSeqnameConvention(seq1, seq2))
            }
        }
    }
}

make_subject <- function() {
    new("GRanges",
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1),
        strand = Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
        values = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

make_query <- function() {
    GRangesList(nomatch = GRanges(seqnames = "chr1",
                                  ranges = IRanges(start=5, end=10),
                                  strand = "+"),
                onematch = GRanges(seqnames = "chr3",
                                   ranges = IRanges(start=2, end=7),
                                   strand = "-"),
                twomatch = GRanges(seqnames = "chr1",
                                   ranges = IRanges(start=1, end=5),
                                   strand = "-"))
}

test_findOverlaps_no_overlaps_returns_empty_matches <- function()
{
    query <- make_query()
    subject <- make_subject()
    ranges(subject) <- shift(ranges(subject), 1000L)

    ## multiple = TRUE
    expect <-
      new("RangesMatching",
          matchMatrix = matrix(integer(),  byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 10L))
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, multiple = TRUE, type = type)
        checkEquals(expect, ans)
    }

    ## multiple = FALSE
    expect <- rep(NA_integer_, length(query))
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, multiple = FALSE, type = type)
        checkEquals(expect, ans)
    }
}

test_findOverlaps_empty_query <- function()
{
    query <- new("GRangesList")
    subject <- make_subject()

    ## multiple = TRUE
    expect <-
      new("RangesMatching",
          matchMatrix = matrix(integer(), byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
            DIM = c(0L, 10L))
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, multiple = TRUE, type = type)
        checkEquals(expect, ans)
    }

    ## multiple = FALSE
    expect <- integer()
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, multiple = FALSE, type = type)
        checkEquals(expect, ans)
    }
}

test_findOverlaps_empty_subject <- function()
{
    query <- make_query()
    subject <- new("GRanges")

    ## multiple = TRUE
    expect <-
      new("RangesMatching",
          matchMatrix = matrix(integer(), byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 0L))
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, multiple = TRUE, type = type)
        checkEquals(expect, ans)
    }

    ## multiple = FALSE
    expect <- rep(NA_integer_, length(query))
    for (type in c("any", "start", "end")) {
        ans <- findOverlaps(query, subject, multiple = FALSE, type = type)
        checkEquals(expect, ans)
    }
}

test_findOverlaps_zero_one_two_matches <- function()
{
    query <- make_query()
    subject <- make_subject()

    ## multiple = TRUE
    expectAny <- 
      new("RangesMatching",
          matchMatrix = matrix(c(2L, 7L, 3L, 1L, 3L, 5L),
                               byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 10L))
    expectStart <- 
      new("RangesMatching",
          matchMatrix = matrix(c(3L, 1L), byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 10L))
    expectEnd <- 
      new("RangesMatching",
          matchMatrix = matrix(integer(), byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 10L))
    ansAny <- findOverlaps(query, subject, multiple = TRUE, type = "any")
    ansStart <- findOverlaps(query, subject, multiple = TRUE, type = "start")
    ansEnd <- findOverlaps(query, subject, multiple = TRUE, type = "end")
    checkEquals(expectAny, ansAny)
    checkEquals(expectStart, ansStart)
    checkEquals(expectEnd, ansEnd)

    ## multiple = FALSE
    expectAny <- c(NA_integer_, 7L, 1L)
    expectStart <- c(NA_integer_, NA_integer_, 1L)
    expectEnd <- c(NA_integer_, NA_integer_, NA_integer_)
    ansAny <- findOverlaps(query, subject, multiple = FALSE, type = "any")
    ansStart <- findOverlaps(query, subject, multiple = FALSE, type = "start")
    ansEnd <- findOverlaps(query, subject, multiple = FALSE, type = "end")
    checkEquals(expectAny, ansAny)
    checkEquals(expectStart, ansStart)
    checkEquals(expectEnd, ansEnd)
}

test_findOverlaps_multimatch_within_one_query <- function()
{
    query <- make_query()
    query[[3L]] <- c(query[[3L]], query[[3L]])
    subject <- make_subject()

    ## multiple = TRUE
    expectAny <- 
      new("RangesMatching",
          matchMatrix = matrix(c(2L, 7L, 3L, 1L, 3L, 5L),
                               byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 10L))
    expectStart <- 
      new("RangesMatching",
          matchMatrix = matrix(c(3L, 1L), byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 10L))
    expectEnd <- 
      new("RangesMatching",
          matchMatrix = matrix(integer(), byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 10L))
    ansAny <- findOverlaps(query, subject, multiple = TRUE, type = "any")
    ansStart <- findOverlaps(query, subject, multiple = TRUE, type = "start")
    ansEnd <- findOverlaps(query, subject, multiple = TRUE, type = "end")
    checkEquals(expectAny, ansAny)
    checkEquals(expectStart, ansStart)
    checkEquals(expectEnd, ansEnd)

    ## multiple = FALSE
    expectAny <- c(NA_integer_, 7L, 1L)
    expectStart <- c(NA_integer_, NA_integer_, 1L)
    expectEnd <- c(NA_integer_, NA_integer_, NA_integer_)
    ansAny <- findOverlaps(query, subject, multiple = FALSE, type = "any")
    ansStart <- findOverlaps(query, subject, multiple = FALSE, type = "start")
    ansEnd <- findOverlaps(query, subject, multiple = FALSE, type = "end")
    checkEquals(expectAny, ansAny)
    checkEquals(expectStart, ansStart)
    checkEquals(expectEnd, ansEnd)
}

test_findOverlaps_missing_strand <- function()
{
    query <- make_query()
    subject <- make_subject()

    query@unlistData@strand <- Rle(strand(c("*", NA, "-")))

    ## multiple = TRUE
    expectAny <- 
      new("RangesMatching",
          matchMatrix = matrix(c(1L, 1L, 1L, 5L, 1L, 6L, 2L, 7L, 3L, 1L, 3L, 5L),
                               byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 10L))
    expectStart <- 
      new("RangesMatching",
          matchMatrix = matrix(c(1L, 5L, 3L, 1L), byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 10L))
    expectEnd <- 
      new("RangesMatching",
          matchMatrix = matrix(c(1L, 1L, 1L, 5L, 1L, 6L),
                               byrow = TRUE, ncol = 2L,
                               dimnames = list(NULL, c("query", "subject"))),
          DIM = c(3L, 10L))
    ansAny <- findOverlaps(query, subject, multiple = TRUE, type = "any")
    ansStart <- findOverlaps(query, subject, multiple = TRUE, type = "start")
    ansEnd <- findOverlaps(query, subject, multiple = TRUE, type = "end")
    checkEquals(expectAny, ansAny)
    checkEquals(expectStart, ansStart)
    checkEquals(expectEnd, ansEnd)

    # multiple = FALSE
    expectAny <- c(1L, 7L, 1L)
    expectStart <- c(5L, NA_integer_, 1L)
    expectEnd <- c(1L, NA_integer_, NA_integer_)
    ansAny <- findOverlaps(query, subject, multiple = FALSE, type = "any")
    ansStart <- findOverlaps(query, subject, multiple = FALSE, type = "start")
    ansEnd <- findOverlaps(query, subject, multiple = FALSE, type = "end")
    checkEquals(expectAny, ansAny)
    checkEquals(expectStart, ansStart)
    checkEquals(expectEnd, ansEnd)
}
