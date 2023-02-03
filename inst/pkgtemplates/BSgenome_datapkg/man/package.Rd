\name{@PKGNAME@}
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
@SRCDATAFILES@
  }

  See \code{?\link[BSgenome]{BSgenomeForge}} and the BSgenomeForge
  vignette (\code{vignette("BSgenomeForge")}) in the \pkg{BSgenome}
  software package for how to create a BSgenome data package.
}

\author{@AUTHOR@}

\seealso{
  \itemize{
    \item \link[BSgenome]{BSgenome} objects in the \pkg{BSgenome}
          software package.

    \item The \code{\link[GenomeInfoDb]{seqinfo}} getter and
          \link[GenomeInfoDb]{Seqinfo} objects in the \pkg{GenomeInfoDb}
          package.

    \item The \code{\link[GenomeInfoDb]{seqlevelsStyle}} getter and
          setter in the \pkg{GenomeInfoDb} package.

    \item \link[Biostrings]{DNAString} objects in the \pkg{Biostrings}
          package.

    \item The \code{\link[BSgenome]{available.genomes}} function
          in the \pkg{BSgenome} software package.

    \item The BSgenomeForge vignette (\code{vignette("BSgenomeForge")})
          in the \pkg{BSgenome} software package for how to create a BSgenome
          data package.
  }
}

\examples{
@PKGNAME@
bsg <- @PKGNAME@
head(seqlengths(bsg))
seqinfo(bsg)

@PKGEXAMPLES@

## ---------------------------------------------------------------------
## Genome-wide motif searching
## ---------------------------------------------------------------------

## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
