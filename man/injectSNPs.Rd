\name{injectSNPs}

\alias{injectSNPs}
\alias{available.SNPs}


\title{injectSNPs}

\description{
   Inject SNPs from a package of SNP locations into a genome.
}

\usage{
  injectSNPs(bsgenome, SNPlocs_pkgname)
  available.SNPs()
}

\arguments{
  \item{bsgenome}{
    [TODO]
  }
  \item{SNPlocs_pkgname}{
    [TODO]
  }
}

\details{
  [TODO]
}

\value{
  [TODO]
}

\author{H. Pages}

\seealso{
  \link{BSgenome-class}
}

\examples{
  ## Look at the SNP locations packages currently available:
  available.SNPs()

  if (interactive()) {
    ## Make your choice and install with:
    source("http://bioconductor.org/biocLite.R")
    biocLite("SNPlocs.Hsapiens.dbSNP.20071016")

    ## Inject into the Human genome:
    library(BSgenome.Hsapiens.UCSC.hg18)
    Hsapiens
    HsWithSNPs <- injectSNPs(Hsapiens, "SNPlocs.Hsapiens.dbSNP.20071016")
    HsWithSNPs
    alphabetFrequency(Hsapiens$chr1)
    alphabetFrequency(HsWithSNPs$chr1)
  }
}

\keyword{manip}
