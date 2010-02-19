## TODO: support strand "*", NA
## TODO: add tests for more complex overlaps
## TODO: performance evaluation and improvement

findOverlaps <- BSgenome:::findOverlaps1

make_target <- function() {
    new("GRanges",
        seqnames = Rle(c("chr1", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
        ranges = IRanges(1:10, width = 10:1),
        strand = Rle(strand(c("-", "+", "+", "-", "-", "-")), c(1, 2, 1, 1, 3, 2)),
        values = DataFrame(score = 1:10, GC = seq(1, 0, length=10)))
}

make_pattern <- function() {
    GRangesList(
                nomatch = GRanges(seqnames = "chr1",
                ranges = IRanges(start=5, end=10), strand = "+"),

                onematch = GRanges( seqnames = "chr3",
                ranges = IRanges(start=2, end=7), strand = "-"),

                twomatch = GRanges(seqnames = "chr1",
                ranges = IRanges(start=1, end=5), strand = "-"))
}

test_findOverlaps_small <- function()
{
    pattern <- make_pattern()
    target <- make_target()
    ans <- findOverlaps(pattern, target)

    expect <- new("RangesMatching",
                  matchMatrix = matrix(c(2L, 7L, 3L, 1L, 3L, 5L),
                  byrow = TRUE, ncol = 2L,
                  dimnames = list(NULL, c("query", "subject"))),
                  DIM = c(10L, 3L))

    checkEquals(expect, ans)
}
