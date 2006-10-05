\name{BSgenome-class}
\docType{class}

% Classes
\alias{class:BSgenome}
\alias{BSgenome-class}

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

\section{Standard generic methods}{
  In the code snippets below,
  \code{x} is a BSgenome object
  and \code{name} is the name of a chromosome (character-string).

  \describe{
    \item{}{
      \code{length(x)}:
      Used to get the length of a BSgenome, i.e., the number
      of \link[Biostrings]{BStringViews} objects that it contains,
      i.e., the number of FASTA files that were used to create this
      BSgenome object (each \link[Biostrings]{BStringViews} object
      actually corresponds to one file).
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
