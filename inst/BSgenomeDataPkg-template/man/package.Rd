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
  \link[BSgenome]{BSgenomeForge}
}

\examples{
@PKGEXAMPLES@

if ("AGAPS" \%in\% masknames(@BSGENOMEOBJNAME@)) {

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

  for (seqname in seqnames(@BSGENOMEOBJNAME@)) {
    cat("Checking sequence", seqname, "... ")
    seq <- @BSGENOMEOBJNAME@[[seqname]]
    checkOnlyNsInGaps(seq)
    cat("OK\n")
  }
}

## See the GenomeSearching vignette in the BSgenome software
## package for some examples of genome-wide motif searching using
## Biostrings and the BSgenome data packages:
if (interactive())
    vignette("GenomeSearching", package="BSgenome")
}

\keyword{package}
\keyword{data}
