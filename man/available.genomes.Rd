\name{available.genomes}

\alias{available.genomes}

\title{available.genomes}

\description{
   Gets the list of Biostrings-based genome data packages
   (aka the BSgenome data packages) currently available.
}

\usage{
  available.genomes()
}

\details{
  A BSgenome data package contains the full genome for a given organism.
  Its name has 4 parts separated by a dot (e.g. BSgenome.Celegans.UCSC.ce2).
  The 1st part is always BSgenome, the 2nd part is the name of the organism
  (abbreviated), the 3rd part is the name of the organisation who assembled
  the genome and the 4th part is the release string or number used by this
  organisation for this genome.
  The package always contains a \link{BSgenome} object named like the second
  part of the package name (e.g. Celegans in the case of BSgenome.Celegans.UCSC.ce2)
  where all the sequences for this genome are stored.
}

\value{
  A character vector containing the names of the BSgenome data packages
  currently available.  
}

\author{H. Pages}

\seealso{
  \link{BSgenome-class}
}

\examples{
  # Look at the BSgenome data packages currently available:
  available.genomes()

  # Make your choice and install with:
  source("http://bioconductor.org/biocLite.R")
  biocLite("BSgenome.Celegans.UCSC.ce2")

  # Have a coffee ;-)

  # Load the package and display the index of sequences for this genome:
  library(BSgenome.Celegans.UCSC.ce2)
  Celegans
}

\keyword{manip}

