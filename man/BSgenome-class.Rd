\name{BSgenome-class}
\docType{class}

% Classes
\alias{class:BSgenome}
\alias{BSgenome-class}

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
      Each single sequence is stored in a \link[Biostrings]{DNAString}
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
  \link[Biostrings]{matchPattern}, \link[Biostrings]{BStringViews},
  \link[Biostrings]{DNAString}
}

\examples{
  library(BSgenome.Celegans.UCSC.ce2)
  length(Celegans)
  Celegans

  ## Access chromosome V (BStringViews object) with
  Celegans[["chrV"]]
  ## or with
  Celegans$chrV

  ## Note that the chromosome V data are not loaded into memory until the
  ## first time you try to access them (what gets loaded is the serialized
  ## BStringViews instance stored in the "chrV.rda" file)
  ## To unload chromosome V (if you need to free some memory)
  unload(Celegans, "chrV")

  ## Note that assignment of BStringViews objects does NOT copy the data
  c <- Celegans[["chrV"]]

  ## But then, clean unloading is a bit more complicated
  unload(Celegans, "chrV")
  remove(c)
  ## More generally, you need to keep track of all BStringViews objects
  ## referencing the same data if you want to be able to unload these data
  ## from memory
}

\keyword{methods}
\keyword{classes}
