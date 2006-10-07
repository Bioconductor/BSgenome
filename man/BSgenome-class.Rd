\name{BSgenome-class}
\docType{class}

% Classes
\alias{class:BSgenome}
\alias{BSgenome-class}

% Constructor-like functions and generics:
\alias{initialize,BSgenome-method}

% Standard generic methods:
\alias{names,BSgenome-method}
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

\section{Standard generic methods}{
  In the code snippets below,
  \code{x} is a BSgenome object
  and \code{name} is the name of a chromosome (character-string).

  \describe{
    \item{}{
      \code{names(x)}:
      Returns the index of \link[Biostrings]{BStringViews} objects
      contained in \code{x}.
      Each \link[Biostrings]{BStringViews} object in \code{x} was created
      from a single FASTA file. The names returned by \code{names(x)}
      are usually the names of those source files but eventually
      a common prefix or suffix was removed in order to keep them as
      short as possible.
    }
    \item{}{
      \code{length(x)}:
      Returns the length of \code{x}, i.e., the number
      of \link[Biostrings]{BStringViews} objects that it contains.
      This is the same as \code{length(names(x))}.
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
  and \code{name} is the name of a chromosome (character-string).

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
