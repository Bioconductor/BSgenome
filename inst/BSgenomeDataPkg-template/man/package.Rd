\name{@BSGENOMEOBJNAME@}
\docType{package}

\alias{@PKGNAME@-package}
\alias{@PKGNAME@}
\alias{@BSGENOMEOBJNAME@}


\title{@PKGTITLE@}

\description{
  @PKGDESCRIPTION@
}

\details{
  @PKGDETAILS@
}

\note{
  This BSgenome data package was made from the following source data files:
  \preformatted{
@SRCDATAFILES1@
@SRCDATAFILES2@
  }
  See \code{?\link[BSgenome]{BSgenomeForge}} and the BSgenomeForge
  vignette (\code{vignette("BSgenomeForge")}) in the BSgenome software
  package for how to make a BSgenome data package.
}

\author{@AUTHOR@}

\seealso{
  \link[BSgenome]{BSgenome-class},
  \link[Biostrings]{DNAString-class},
  \code{\link[BSgenome]{available.genomes}},
  \link[BSgenome]{BSgenomeForge},
  \link{@PKGNAME@-datasets}
}

\examples{
@PKGEXAMPLES@

## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and a BSgenome data package:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
