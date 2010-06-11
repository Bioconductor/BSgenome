\name{getSeq}

\alias{getSeq}
\alias{getSeq,BSgenome-method}

\title{getSeq}

\description{
   A function for extracting a set of sequences (or subsequences) from
   a \link{BSgenome} object or other sequence container. This man page 
   specifically documents the method for \link{BSgenome} objects.
}

\usage{
getSeq(x, ...)

\S4method{getSeq}{BSgenome}(x, names, start=NA, end=NA, width=NA,
                 strand="+", as.character=TRUE) 
}

\arguments{
  \item{x}{
    A \link{BSgenome} object.
    See the \code{\link{available.genomes}} function for how
    to install a genome.
  }
  \item{names}{
    The names of the sequences in \code{x} from where to extract the
    subsequences, or a \link[GenomicRanges]{GRanges} object,
    or a \link[IRanges]{RangedData} object,
    or a named \link[IRanges]{RangesList} object,
    or a named \link[IRanges]{Ranges} object.
    The \link[IRanges]{RangesList} or \link[IRanges]{Ranges} object
    must be named according to the sequences in \code{x} from where
    to extract the subsequences.

    If \code{names} is missing, then \code{seqnames(x)} is used.
    
    See \code{?`\link{seqnames,BSgenome-method}`} and
    \code{?`\link{mseqnames,BSgenome-method}`} to get 
    the list of single sequences and multiple sequences (respectively)
    contained in \code{x}.
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
  then it is also recycled to length L.

  Here is how the lookup between the names passed to the \code{names}
  argument and the sequences in \code{x} is performed.
  For each \code{name} in \code{names}:
  \itemize{
  \item (1): If \code{x} contains a single sequence with that name
        then this sequence is used for extraction;
  \item (2): Otherwise the names of all the elements in all the
        multiple sequences are searched: \code{name} is treated
        as a regular expression and \code{\link[base]{grep}} is
        used for this search. If exactly one sequence is found,
        then it is used for extraction, otherwise an error is raised.
  }
}

\value{
  A character vector of length L when \code{as.character=TRUE}.

  A \link[Biostrings]{DNAString} or \link[Biostrings]{DNAStringSet}
  object when \code{as.character=FALSE}.
  More precisely the returned value is a \link[Biostrings]{DNAString}
  object if L = 1 and \code{names} is not a \link[GenomicRanges]{GRanges},
  \link[IRanges]{RangedData}, \link[IRanges]{RangesList} or
  \link[IRanges]{Ranges} object.
  Otherwise it's a \link[Biostrings]{DNAStringSet} object.
}

\note{
  The default value for the \code{as.character} argument is \code{TRUE}
  even though \code{getSeq} is much more efficient (in terms of speed
  and memory usage) when used with \code{as.character=FALSE}.

  Be aware that using \code{as.character=TRUE} can be very inefficient
  when extracting a high number of sequences (hundreds of thousands)
  or long sequences (> 1 million letters).

  For this reason the plan is to change the default value to \code{FALSE}
  in the near future.

  Note that the masks in \code{x}, if any, are always ignored. In other
  words, masked regions in the genome are extracted in the same way as
  unmasked regions (this is achieved by dropping the masks before extraction).
  See \code{?`\link[Biostrings]{MaskedXString-class}`} for more
  information about masked sequences.
}

\author{H. Pages; improvements suggested by Matt Settles and others}

\seealso{
  \code{\link{available.genomes}},
  \link{BSgenome-class},
  \link{seqnames,BSgenome-method},
  \link{mseqnames,BSgenome-method},
  \link{[[,BSgenome-method},
  \link[Biostrings]{DNAString-class},
  \link[Biostrings]{MaskedDNAString-class},
  \link[GenomicRanges]{GRanges-class},
  \link[IRanges]{RangesList-class},
  \link[IRanges]{Ranges-class},
  \link[IRanges]{subseq,DataTable},
  \link[IRanges]{subseq,Sequence},
  \code{\link[base]{grep}}
}

\examples{
  ## ---------------------------------------------------------------------
  ## A. EXTRACTING A SMALL NUMBER OF SEQUENCES FROM THE C. elegans GENOME
  ## ---------------------------------------------------------------------

  # Load the Caenorhabditis elegans genome (UCSC Release ce2):
  library(BSgenome.Celegans.UCSC.ce2)

  # Look at the index of sequences:
  Celegans

  # Get chromosome V as a DNAString object:
  getSeq(Celegans, "chrV", as.character=FALSE)
  # which is in fact the same as doing:
  Celegans$chrV

  # Never try this:
  #getSeq(Celegans, "chrV")
  # or this (even worse):
  #getSeq(Celegans)

  # Get the first 20 bases of each chromosome:
  getSeq(Celegans, end=20)

  # Get the last 20 bases of each chromosome:
  getSeq(Celegans, start=-20)

  # Extracting small sequences from different chromosomes:
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

  # Get the "NM_058280_up_1000" sequence (belongs to the upstream1000
  # multiple sequence) as a character string:
  s1 <- getSeq(Celegans, "NM_058280_up_1000")
  # or a DNAString object (more efficient):
  s2 <- getSeq(Celegans, "NM_058280_up_1000", as.character=FALSE)

  getSeq(Celegans, "NM_058280_up_5000", start=-1000) == s1  # TRUE

  getSeq(Celegans, "NM_058280_up_5000",
         start=-1000, as.character=FALSE) == s2  # TRUE

  ## ---------------------------------------------------------------------
  ## B. RANDOMLY EXTRACTING A HIGH NUMBER OF 40-mers FROM THE C. elegans
  ##    GENOME
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
      getSeq(x, gr, as.character=FALSE)
  }
  library(BSgenome.Celegans.UCSC.ce2)
  ## With a density of 1 read every 100 genome bases, the total number of
  ## extracted reads is about 1 million:
  rndreads <- extractRandomReads(Celegans, 0.01, 40)
}

\keyword{manip}
