\name{@PKGNAME@}
\docType{package}

\alias{@PKGNAME@-package}
\alias{@PKGNAME@}

\title{@PKGTITLE@}

\description{
  @PKGDESCRIPTION@
}

\details{
  @PKGDETAILS@
}

\note{
  The masks in this BSgenome data package were made from the following
  source data files:
  \preformatted{
@SRCDATAFILES@
  }

  See \code{?\link[@REFPKGNAME@]{@REFPKGNAME@}} in the
  \pkg{@REFPKGNAME@} package for information about how the sequences
  were obtained.

  See \code{?\link[BSgenome]{BSgenomeForge}} and the BSgenomeForge
  vignette (\code{vignette("BSgenomeForge")}) in the \pkg{BSgenome}
  software package for how to create a BSgenome data package.
}

\author{@AUTHOR@}

\seealso{
  \itemize{
    \item \link[@REFPKGNAME@]{@REFPKGNAME@} in the \pkg{@REFPKGNAME@} package
          for information about how the sequences were obtained.

    \item \link[BSgenome]{BSgenome} objects in the \pkg{BSgenome}
          software package.

    \item The \code{\link[GenomeInfoDb]{seqinfo}} getter and
          \link[GenomeInfoDb]{Seqinfo} objects in the \pkg{GenomeInfoDb}
          package.

    \item The \code{\link[GenomeInfoDb]{seqlevelsStyle}} getter and
          setter in the \pkg{GenomeInfoDb} package.

    \item \link[Biostrings]{MaskedDNAString} objects in the \pkg{Biostrings}
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
mbsg <- @PKGNAME@
head(seqlengths(mbsg))
seqinfo(mbsg)

@PKGEXAMPLES@

if ("AGAPS" \%in\% masknames(mbsg)) {

  ## Check that the assembly gaps contain only Ns:
  checkOnlyNsInGaps <- function(seq)
  {
    ## Replace all masks by the inverted AGAPS mask
    masks(seq) <- gaps(masks(seq)["AGAPS"])
    unique_letters <- uniqueLetters(seq)
    if (any(unique_letters != "N"))
        stop("assembly gaps contain more than just Ns")
  }

  ## A message will be printed each time a sequence is removed
  ## from the cache:
  options(verbose=TRUE)

  for (seqname in seqnames(mbsg)) {
    cat("Checking sequence", seqname, "... ")
    seq <- mbsg[[seqname]]
    checkOnlyNsInGaps(seq)
    cat("OK\n")
  }
}

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
