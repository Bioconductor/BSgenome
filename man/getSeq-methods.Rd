\name{getSeq-methods}

\alias{getSeq,BSgenome-method}

\title{getSeq method for BSgenome objects}

\description{
   A \code{\link[Biostrings]{getSeq}} method for extracting a set of
   sequences (or subsequences) from a \link{BSgenome} object.
}

\usage{
\S4method{getSeq}{BSgenome}(x, names, start=NA, end=NA, width=NA,
                 strand="+", as.character=FALSE) 
}

\arguments{
  \item{x}{
    A \link{BSgenome} object.
    See the \code{\link{available.genomes}} function for how
    to install a genome.
  }
  \item{names}{
    A character vector containing the names of the sequences
    in \code{x} where to get the subsequences from,
    or a \link[GenomicRanges]{GRanges} object,
    or a \link[IRanges]{RangedData} object,
    or a named \link[IRanges]{RangesList} object,
    or a named \link[IRanges]{Ranges} object.
    The \link[IRanges]{RangesList} or \link[IRanges]{Ranges} object
    must be named according to the sequences
    in \code{x} where to get the subsequences from.

    If \code{names} is missing, then \code{seqnames(x)} is used.
    
    See \code{?`\link{BSgenome-class}`} for details on how to get
    the lists of single sequences and multiple sequences (respectively)
    contained in a \link{BSgenome} object.
  }
  \item{start, end, width}{
    Vector of integers (eventually with NAs) specifying the locations
    of the subsequences to extract.
    These are not needed (and therefore it's an error to supply them)
    when \code{names} is a \link[GenomicRanges]{GRanges},
    \link[IRanges]{RangedData}, \link[IRanges]{RangesList}
    or \link[IRanges]{Ranges} object.
  }
  \item{strand}{
    A vector containing \code{"+"}s or/and \code{"-"}s.
    This is not needed (and therefore it's an error to supply it)
    when \code{names} is a \link[GenomicRanges]{GRanges} object
    or a \link[IRanges]{RangedData} object with a strand column.
  }
  \item{as.character}{
    \code{TRUE} or \code{FALSE}. Should the extracted sequences
    be returned in a standard character vector?
  }
  \item{...}{Additional arguments. (Currently ignored.)}
}

\details{
  L, the number of sequences to extract, is determined as follow:
  \itemize{
  \item If \code{names} is a \link[GenomicRanges]{GRanges} or
        \link[IRanges]{Ranges} object then L = \code{length(names)}.
  \item If \code{names} is a \link[IRanges]{RangedData} object
        then L = \code{nrow(names)}.
  \item If \code{names} is a \link[IRanges]{RangesList} object
        then L = \code{length(unlist(names))}.
  \item Otherwise, L is the length of the longest of \code{names},
        \code{start}, \code{end} and \code{width} and all these
        arguments are recycled to this length.
        \code{NA}s and negative values in these 3 arguments are
        solved according to the rules of the SEW (Start/End/Width)
        interface (see \code{?\link[IRanges]{solveUserSEW}} for
        the details).
  }

  If \code{names} is neither a \link[GenomicRanges]{GRanges} object
  or a \link[IRanges]{RangedData} object with a strand column,
  then the \code{strand} argument is also recycled to length L.

  Here is how the lookup between the names passed to the \code{names}
  argument and the sequences in \code{x} is performed.
  For each \code{name} in \code{names}:
  \itemize{
  \item (1): If \code{x} contains a single sequence with that name
        then this sequence is used for extraction;
  \item (2): Otherwise the names of all the elements in all the
        multiple sequences are searched. If the \code{names} argument
        is a character vector then \code{name} is treated as a regular
        expression and \code{\link[base]{grep}} is used for this search,
        otherwise (i.e. when the names are supplied via a higher level
        object like \link[GenomicRanges]{GRanges}) \code{name} must match
        exactly the name of the sequence. If exactly one sequence is found,
        then it is used for extraction, otherwise an error is raised.
  }
}

\value{
  A character vector of length L when \code{as.character=TRUE}.

  A \link[Biostrings]{DNAString} or \link[Biostrings]{DNAStringSet}
  object when \code{as.character=FALSE} (the default).
  More precisely the returned value is a \link[Biostrings]{DNAString}
  object if L = 1 and \code{names} is not a \link[GenomicRanges]{GRanges},
  \link[IRanges]{RangedData}, \link[IRanges]{RangesList} or
  \link[IRanges]{Ranges} object.
  Otherwise it's a \link[Biostrings]{DNAStringSet} object.
}

\note{
  Be aware that using \code{as.character=TRUE} can be very inefficient
  when extracting a "big" amount of DNA sequences (e.g. millions of
  short sequences or a small number of very long sequences).

  Note that the masks in \code{x}, if any, are always ignored. In other
  words, masked regions in the genome are extracted in the same way as
  unmasked regions (this is achieved by dropping the masks before extraction).
  See \code{?`\link[Biostrings]{MaskedDNAString-class}`} for more
  information about masked DNA sequences.
}

\author{H. Pages; improvements suggested by Matt Settles and others}

\seealso{
  \code{\link[Biostrings]{getSeq}},
  \code{\link{available.genomes}},
  \link{BSgenome-class},
  \link[Biostrings]{DNAString-class},
  \link[Biostrings]{DNAStringSet-class},
  \link[Biostrings]{MaskedDNAString-class},
  \link[GenomicRanges]{GRanges-class},
  \link[IRanges]{RangedData-class},
  \link[IRanges]{RangesList-class},
  \link[IRanges]{Ranges-class},
  \code{\link[base]{grep}}
}

\examples{
## ---------------------------------------------------------------------
## A. SIMPLE EXAMPLES
## ---------------------------------------------------------------------

## Load the Caenorhabditis elegans genome (UCSC Release ce2):
library(BSgenome.Celegans.UCSC.ce2)

## Look at the index of sequences:
Celegans

## Get chromosome V as a DNAString object:
getSeq(Celegans, "chrV")
## which is in fact the same as doing:
Celegans$chrV

\dontrun{
  ## Never try this:
  getSeq(Celegans, "chrV", as.character=TRUE)
  ## or this (even worse):
  getSeq(Celegans, as.character=TRUE)
}

## Get the first 20 bases of each chromosome:
getSeq(Celegans, end=20)

## Get the last 20 bases of each chromosome:
getSeq(Celegans, start=-20)

## Get the "NM_058280_up_1000" sequence (belongs to the upstream1000
## multiple sequence) as a DNAString object:
s1 <- getSeq(Celegans, "NM_058280_up_1000")
stopifnot(identical(getSeq(Celegans, "NM_058280_up_5000", start=-1000), s1))

\dontrun{
  ## Fails because there is more than one sequence across
  ## Celegans$upstream1000, Celegans$upstream2000 and Celegans$upstream5000
  ## with "NM_058280" in its name:
  getSeq(Celegans, "NM_058280")

  ## Fails because there is no sequence named exactly "NM_058280":
  getSeq(Celegans, "^NM_058280$")
}

## ---------------------------------------------------------------------
## B. EXTRACTING SMALL SEQUENCES FROM DIFFERENT CHROMOSOMES
## ---------------------------------------------------------------------

myseqs <- data.frame(
  chr=c("chrI", "chrX", "chrM", "chrM", "chrX", "chrI", "chrM", "chrI"),
  start=c(NA, -40, 8510, 301, 30001, 9220500, -2804, -30),
  end=c(50, NA, 8522, 324, 30011, 9220555, -2801, -11),
  strand=c("+", "-", "+", "+", "-", "-", "+", "-")
)
getSeq(Celegans, myseqs$chr,
       start=myseqs$start, end=myseqs$end)
getSeq(Celegans, myseqs$chr,
       start=myseqs$start, end=myseqs$end, strand=myseqs$strand)

## ---------------------------------------------------------------------
## C. USING A GRanges OBJECT
## ---------------------------------------------------------------------

gr1 <- GRanges(seqnames=c("chrI", "chrI", "chrM"),
               ranges=IRanges(start=101:103, width=9))
gr1  # all strand values are "*"
getSeq(Celegans, gr1)  # treats strand values as if they were "+"

strand(gr1)[] <- "-"
getSeq(Celegans, gr1)

strand(gr1)[1] <- "+"
getSeq(Celegans, gr1)

strand(gr1)[2] <- "*"
if (interactive())
  getSeq(Celegans, gr1)  # Error: cannot mix "*" with other strand values

gr2 <- GRanges(seqnames=c("chrM", "NM_058280_up_1000"),
               ranges=IRanges(start=103:102, width=9))
gr2
if (interactive()) {
  ## Because the sequence names are supplied via a GRanges object, they
  ## are not treated as regular expressions:
  getSeq(Celegans, gr2)  # Error: sequence NM_058280_up_1000 not found
}

## ---------------------------------------------------------------------
## D. EXTRACTING A HIGH NUMBER OF RANDOM 40-MERS FROM A GENOME
## ---------------------------------------------------------------------

extractRandomReads <- function(x, density, readlength)
{
    if (!is.integer(readlength))
        readlength <- as.integer(readlength)
    start <- lapply(seqnames(x),
                    function(name)
                    {
                      seqlength <- seqlengths(x)[name]
                      sample(seqlength - readlength + 1L,
                             seqlength * density,
                             replace=TRUE)
                    })
    names <- rep.int(seqnames(x), elementLengths(start))
    ranges <- IRanges(start=unlist(start), width=readlength)
    strand <- strand(sample(c("+", "-"), length(names), replace=TRUE))
    gr <- GRanges(seqnames=names, ranges=ranges, strand=strand)
    getSeq(x, gr)
}

## With a density of 1 read every 100 genome bases, the total number of
## extracted 40-mers is about 1 million:
rndreads <- extractRandomReads(Celegans, 0.01, 40)

## Notes:
## - The short sequences in 'rndreads' can be seen as the result of a
##   simulated high-throughput sequencing experiment. A non-realistic
##   one though because:
##     (a) It assumes that the underlying technology is perfect (the
##         generated reads have no technology induced errors).
##     (b) It assumes that the sequenced genome is exactly the same as
##         the reference genome.
##     (c) The simulated reads can contain IUPAC ambiguity letters only
##         because the reference genome contains them. In a real
##         high-throughput sequencing experiment, the sequenced genome
##         of course doesn't contain those letters, but the sequencer
##         can introduce them in the generated reads to indicate
##         ambiguous base-calling.
## - Those reads are coming from the plus and minus strands of the
##   chromosomes.
## - With a density of 0.01 and the reads being only 40-base long, the
##   average coverage of the genome is only 0.4 which is low. The total
##   number of reads is about 1 million and it takes less than 10 sec.
##   to generate them.
## - A higher coverage can be achieved by using a higher density and/or
##   longer reads. For example, with a density of 0.1 and 100-base reads
##   the average coverage is 10. The total number of reads is about 10
##   millions and it takes less than 1 minute to generate them.
## - Those reads could easily be mapped back to the reference by using
##   an efficient matching tool like matchPDict() for performing exact
##   matching (see ?matchPDict for more information). Typically, a
##   small percentage of the reads (4 to 5\% in our case) will hit the
##   reference at multiple locations. This is especially true for such
##   short reads, and, in a lower proportion, is still true for longer
##   reads, even for reads as long as 300 bases.

## ---------------------------------------------------------------------
## E. SEE THE BSgenome CACHE IN ACTION
## ---------------------------------------------------------------------

options(verbose=TRUE)
first20 <- getSeq(Celegans, end=20)
first20
gc()
stopifnot(length(ls(Celegans@.seqs_cache)) == 0L)
## One more gc() call is needed in order to see the amount of memory in
## used after all the chromosomes have been removed from the cache:
gc()
}

\keyword{manip}
