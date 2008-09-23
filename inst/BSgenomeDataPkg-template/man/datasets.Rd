\name{@PKGNAME@-datasets}
\docType{data}

\alias{@PKGNAME@-datasets}
\alias{seqlengths}
@MASKSDATASETALIASES@

\title{@PKGNAME@ datasets}

\description{
  Data sets included in the @PKGNAME@ package in addition to the
  sequence data (and related to the sequence data).
}

\usage{
## Load the lengths of the single sequences contained in the package:
data(seqlengths)

@MASKSDATASETUSAGES@
}

\author{@AUTHOR@}

\seealso{
  \link{@PKGNAME@},
  \link[BSgenome]{BSgenomeForge}
}

\examples{
  ## You need to specify the 'package' argument if @PKGNAME@
  ## is not the first package in your search path:
  data(seqlengths, package="@PKGNAME@")
  seqlengths  # same as seqlengths(@BSGENOMEOBJNAME@)
}

\keyword{package}
\keyword{data}
\keyword{datasets}
