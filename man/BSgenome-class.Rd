\name{BSgenome-class}
\docType{class}

% Classes
\alias{class:BSgenome}
\alias{BSgenome-class}
\alias{BSgenome}

% Accesor methods:
\alias{seqnames}
\alias{seqnames,BSgenome-method}
\alias{mseqnames}
\alias{mseqnames,BSgenome-method}
\alias{names,BSgenome-method}

% Constructor-like functions and generics:
\alias{initialize,BSgenome-method}

% Standard generic methods:
\alias{show,BSgenome-method}
\alias{length,BSgenome-method}
\alias{[[,BSgenome-method}
\alias{[[<-,BSgenome-method}
\alias{$,BSgenome-method}

% Other functions and generics:
\alias{unload}
\alias{unload,BSgenome-method}


\title{The BSgenome class}

\description{
  A container for the complete genome sequence of a given specie.
}

\details{
    [TODO: Put some details here]
}

\section{Accesor methods}{
  In the code snippets below,
  \code{x} is a BSgenome object.

  \describe{
    \item{}{
      \code{seqnames(x)}:
      Returns the index of the single sequences contained in \code{x}.
      Each single sequence is stored in a \link[Biostrings]{BString} (or derived)
      object and comes from a source file (FASTA) with a single record.
      The names returned by \code{seqnames(x)} usually reflect the names
      of those source files but a common prefix or suffix was eventually
      removed in order to keep them as short as possible.
    }
    \item{}{
      \code{mseqnames(x)}:
      Returns the index of the multiple sequences contained in \code{x}.
      Each multiple sequence is stored in a \link[Biostrings]{BStringViews}
      object and comes from a source file (FASTA) with multiple records.
      The names returned by \code{mseqnames(x)} usually reflect the names
      of those source files but a common prefix or suffix was eventually
      removed in order to keep them as short as possible.
    }
    \item{}{
      \code{names(x)}:
      Returns the index of all sequences contained in \code{x}.
      This is the same as \code{c(seqnames(x), mseqnames(x))}.
    }
  }
}

\section{Standard generic methods}{
  In the code snippets below,
  \code{x} is a BSgenome object
  and \code{name} is the name of a sequence (character-string).

  \describe{
    \item{}{
      \code{length(x)}:
      Returns the length of \code{x}, i.e., the number of all sequences
      that it contains. This is the same as \code{length(names(x))}.
    }
    \item{}{
      \code{x[[name]]}:
      [TODO: Document me]
    }
    \item{}{
      \code{x$name}:
      [TODO: Document me]
    }
  }
}

\section{Other functions and generics}{
  In the code snippets below,
  \code{x} is a BSgenome object
  and \code{name} is the name of a sequence (character-string).

  \describe{
    \item{}{
      \code{unload(x, name)}:
      [TODO: Document me]
    }
  }
}

\author{H. Pages}

\seealso{
  \code{\link{available.genomes}},
  \link[Biostrings]{BString},
  \link[Biostrings]{DNAString},
  \link[Biostrings]{BStringViews},
  \code{\link{getSeq}},
  \code{\link[Biostrings]{matchPattern}},
  \code{\link{rm}},
  \code{\link{gc}}
}

\examples{
  library(BSgenome.Celegans.UCSC.ce2)   # This doesn't load the chromosome 
                                        # sequences into memory.
  length(Celegans)                      # Number of sequences in this genome.
  Celegans                              # Displays index of all the sequences
                                        # in this genome.
  mem0 <- gc()["Vcells", "(Mb)"]        # Current amount of data in memory (in
                                        # Mb).
  Celegans[["chrV"]]                    # Loads chromosome V into memory (hence
                                        # takes a long time).
  gc()["Vcells", "(Mb)"] - mem0         # Chromosome V occupies 20Mb of memory.
  Celegans[["chrV"]]                    # Much faster (sequence is already in
                                        # memory, hence it's not loaded again).
  Celegans$chrV                         # Equivalent to Celegans[["chrV"]].
  class(Celegans$chrV)                  # Chromosome V (like any other
                                        # chromosome sequence) is a DNAString
                                        # object.
  nchar(Celegans$chrV)                  # Its has 20922231 letters (nucleotides).
  x <- Celegans$chrV                    # Very fast because a BString object
                                        # doesn't contain the sequence, only a
                                        # pointer to the sequence, hence chrV
                                        # seq is not duplicated in memory. But
                                        # we now have 2 objects pointing to the
                                        # same place in memory.
  y <- substr(x, 10, 100)               # A 3rd object pointing to chrV seq.
  
  ## We must remove all references to chrV seq if we want the 20Mb of memory
  ## used by it to be freed (note that it can be hard to keep track of all the
  ## references to a given sequence).
  ## IMPORTANT: The 1st reference to this seq (Celegans$chrV) should be removed
  ## last. This is achieved with unload(). All other references are removed by
  ## just removing the referencing object.
  rm(x)
  rm(y)
  unload(Celegans, "chrV")
  gc()["Vcells", "(Mb)"]
}

\keyword{methods}
\keyword{classes}

