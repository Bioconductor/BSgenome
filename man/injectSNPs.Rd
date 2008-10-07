\name{injectSNPs}

\alias{injectSNPs}
\alias{injectSNPs,BSgenome-method}
\alias{SNPlocs_pkg}
\alias{SNPlocs_pkg,BSgenome-method}
\alias{SNPcount}
\alias{SNPcount,BSgenome-method}
\alias{SNPlocs}
\alias{SNPlocs,BSgenome-method}

\alias{available.SNPs}


\title{SNP injection}

\description{
   Inject SNPs from a package of SNP locations into a genome.
}

\usage{
  available.SNPs()
  injectSNPs(x, SNPlocs_pkg)

  ## Related utilities
  SNPlocs_pkg(x)
  SNPcount(x)
  SNPlocs(x, seqname)
}

\arguments{
  \item{x}{
    A \link{BSgenome} object.
  }
  \item{SNPlocs_pkg}{
    The name of a package containing SNP locations for the single
    sequences contained in \code{x}.
    This package must have been previously installed by the user,
    \code{injectSNPs} won't try to install it.
  }
  \item{seqname}{
    The name of a single sequence in \code{x}.
  }
}

\value{
  All the functions below have the side effect to load the package
  containing the SNP locations.

  \code{injectSNPs} returns a copy of the original genome \code{x}
  where some or all of the single sequences were altered by injecting
  the SNPs defined in the \code{SNPlocs_pkg} package.

  \code{SNPlocs_pkg}, \code{SNPcount} and \code{SNPlocs} return \code{NULL}
  if no SNPs were injected in \code{x} (i.e. if \code{x} is not a
  \link{BSgenome} object returned by a previous call to \code{injectSNPs}).
  Otherwise \code{SNPlocs_pkg} returns the name of the package from which the
  SNPs were injected, \code{SNPcount} the number of SNPs for each altered
  sequence in \code{x}, and \code{SNPlocs} their locations in the sequence
  whose name is specified by \code{seqname}.
}

\author{H. Pages}

\seealso{
  \link{BSgenome-class}
}

\examples{
  ## Look at the packages of SNP locations currently available:
  available.SNPs()

  if (interactive()) {
    ## Make your choice and install with:
    source("http://bioconductor.org/biocLite.R")
    biocLite("SNPlocs.Hsapiens.dbSNP.20071016")
  }

  ## Inject SNPs from dbSNP into the Human genome:
  library(BSgenome.Hsapiens.UCSC.hg18)
  Hsapiens
  HsWithSNPs <- injectSNPs(Hsapiens, "SNPlocs.Hsapiens.dbSNP.20071016")
  HsWithSNPs
  alphabetFrequency(Hsapiens$chr1)
  alphabetFrequency(HsWithSNPs$chr1)
}

\keyword{manip}
