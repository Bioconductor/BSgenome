\name{BSgenomeForge}

\alias{BSgenomeForge}
\alias{forgeSeqFiles}
\alias{forgeMaskFiles}

\title{The BSgenomeForge functions}

\description{
  The BSgenomeForge functions are intended to be used for making
  "Biostring-based genome data packages" (aka "BSgenome packages")
  like the BSgenome.Hsapiens.UCSC.hg18 package.
  They are not intended to be used by "regular" users of these packages.

  There are currently two BSgenomeForge functions: \code{forgeSeqFiles}
  and \code{forgeMaskFiles}.
}

\usage{
  forgeSeqFiles(srcdir, destdir, names, prefix="", suffix="", comments=NULL, single.seq=TRUE)
  forgeMaskFiles(srcdir, destdir, seqnames, seqdir, nmasks, prefix="", suffix="")
}

\arguments{
  \item{srcdir}{
    [TODO: Document me]
  }
  \item{destdir}{
    [TODO: Document me]
  }
  \item{names}{
    [TODO: Document me]
  }
  \item{prefix}{
    [TODO: Document me]
  }
  \item{suffix}{
    [TODO: Document me]
  }
  \item{comments}{
    [TODO: Document me]
  }
  \item{single.seq}{
    [TODO: Document me]
  }
  \item{seqnames}{
    [TODO: Document me]
  }
  \item{seqdir}{
    [TODO: Document me]
  }
  \item{nmasks}{
    [TODO: Document me]
  }
}

\author{H. Pages}

\keyword{manip}
