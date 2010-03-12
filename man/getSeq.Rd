\name{getSeq}

\alias{getSeq}
\alias{getSeq,BSgenome-method}

\title{getSeq}

\description{
   A convenience function for extracting a set of sequences (or
   subsequences) from a \link[BSgenome:BSgenome-class]{BSgenome} or
   other object. This man page 
   specifically documents the \link[BSgenome:BSgenome-class]{BSgenome} method.
}

\usage{
  getSeq(x, ...)

  \S4method{getSeq}{BSgenome}(x, names, start=NA, end=NA, width=NA,
                strand="+", as.character=TRUE) 
}

\arguments{
  \item{x}{
    A \code{\link{BSgenome}} object.
    See the \code{\link{available.genomes}} function for how
    to install a genome.
  }
  \item{names}{
    The names of the sequences to extract from \code{x}, or a
    \code{RangedData} or \code{RangesList} object.
    If missing, then \code{seqnames(x)} is used.
    
    See \code{?`\link{seqnames,BSgenome-method}`} and
    \code{?`\link{mseqnames,BSgenome-method}`} to get 
    the list of single sequences and multiple sequences (respectively)
    contained in \code{x}.

    Here is how the lookup between the names passed to the \code{names}
    argument and the sequences in \code{x} is performed.
    For each \code{name} in \code{names}:
    (1) if \code{x} contains a single sequence with that name
    then this sequence is returned;
    (2) otherwise the names of all the elements in all the multiple
    sequences are searched: \code{name} is treated as a regular
    expression and \code{\link[base]{grep}} is used for this search.
    If exactly one sequence is found, then it's returned, otherwise an
    error is raised.

    If \code{names} is a \code{RangedData} or \code{RangesList}, the
    \code{space}, \code{start}, and \code{width} are extracted and
    treated as the \code{names}, \code{start} and \code{width}
    arguments, respectively. In the case of a \code{RangedData}, the
    \dQuote{strand} column, if any, is extracted and overrides the
    \code{strand} argument. If there is no \dQuote{strand} column, all
    features are assumed to be on the positive strand. A warning is
    emitted if any of the overridden arguments is passed to the
    function.
  }
  \item{start, end, width}{
    Vector of integers (eventually with NAs). Overridden if \code{names}
    is a \code{RangedData} or \code{RangesList}.
  }
  \item{strand}{
    A vector containing \code{+}s or/and \code{-}s. Overridden if
    \code{names} is a \code{RangedData}.
  }
  \item{as.character}{
    \code{TRUE} or \code{FALSE}. Should the extracted sequences
    be returned in a standard character vector?
  }

  \item{...}{Additional arguments. (Currently ignored.)}
}

\details{
  The \code{names}, \code{start}, \code{end}, \code{width} and \code{strand}
  arguments are expanded cyclically to the length of the longest provided
  none are of zero length.
}

\value{
  A standard character vector when \code{as.character=TRUE}.
  Note that when \code{as.character=TRUE}, then the masks that
  are defined on top of the sequences to extract are ignored
  (i.e. dropped) if any (see \code{?`\link[Biostrings]{MaskedXString-class}`}
  for more information about masked sequences).

  A \link[Biostrings:DNAString-class]{DNAString} or
  \link[Biostrings:MaskedXString-class]{MaskedDNAString} 
  object when \code{as.character=FALSE}.
  Note that \code{as.character=FALSE} is not supported yet when extracting
  more than one sequence.
}

\note{
  Be aware that using \code{as.character=TRUE} can be very inefficient
  when the returned character vector contains very long strings
  (> 1 million letters) or is itself a long vector (> 10000 strings).

  \code{getSeq} is much more efficient when used with
  \code{as.character=FALSE} but this works only for extracting
  one sequence at a time for now.
}

\author{H. Pages; improvements suggested by Matt Settles and others}

\seealso{
  \code{\link{available.genomes}},
  \code{\link{BSgenome-class}},
  \code{\link{seqnames,BSgenome-method}},
  \code{\link{mseqnames,BSgenome-method}},
  \code{\link{[[,BSgenome-method}},
  \code{\link[base]{grep}},
  \code{\link[IRanges:DataTable-class]{subseq,DataTable}},
  \code{\link[IRanges:Sequence-class]{subseq,Sequence}},
  \code{\link[Biostrings:DNAString-class]{DNAString}},
  \code{\link[Biostrings:MaskedXString-class]{MaskedDNAString}}
}

\examples{
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
}

\keyword{manip}
