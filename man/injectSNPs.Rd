\name{injectSNPs}

\alias{class:InjectSNPsHandler}
\alias{InjectSNPsHandler-class}
\alias{InjectSNPsHandler}

\alias{injectSNPs}
\alias{injectSNPs,BSgenome-method}
\alias{SNPlocs_pkgname}
\alias{SNPlocs_pkgname,InjectSNPsHandler-method}
\alias{SNPlocs_pkgname,BSgenome-method}
\alias{SNPcount}
\alias{SNPcount,InjectSNPsHandler-method}
\alias{SNPcount,BSgenome-method}
\alias{SNPlocs}
\alias{SNPlocs,InjectSNPsHandler-method}
\alias{SNPlocs,BSgenome-method}

\alias{available.SNPs}
\alias{installed.SNPs}


\title{SNP injection}

\description{
  Inject SNPs from a SNPlocs data package into a genome.
}

\usage{
injectSNPs(x, SNPlocs_pkgname)

SNPlocs_pkgname(x)
SNPcount(x)
SNPlocs(x, seqname)

## Related utilities
available.SNPs(type=getOption("pkgType"))
installed.SNPs()
}

\arguments{
  \item{x}{
    A \link[BSgenome:BSgenome-class]{BSgenome} object.
  }
  \item{SNPlocs_pkgname}{
    The name of a SNPlocs data package containing SNP information for the
    single sequences contained in \code{x}.
    This package must be already installed (\code{injectSNPs} won't try to
    install it).
  }
  \item{seqname}{
    The name of a single sequence in \code{x}.
  }
  \item{type}{
    Character string indicating the type of package (\code{"source"},
    \code{"mac.binary"} or \code{"win.binary"}) to look for.
  }
}

\value{
  \code{injectSNPs} returns a copy of the original genome \code{x} where some
  or all of the single sequences were altered by injecting the SNPs defined in
  the SNPlocs data package specified thru the \code{SNPlocs_pkgname} argument.
  The SNPs in the altered genome are represented by an IUPAC ambiguity code
  at each SNP location.

  \code{SNPlocs_pkgname}, \code{SNPcount} and \code{SNPlocs} return \code{NULL}
  if no SNPs were injected in \code{x} (i.e. if \code{x} is not a
  \link[BSgenome:BSgenome-class]{BSgenome} object returned by a previous
  call to \code{injectSNPs}). Otherwise \code{SNPlocs_pkgname} returns
  the name of the package from  which the 
  SNPs were injected, \code{SNPcount} the number of SNPs for each altered
  sequence in \code{x}, and \code{SNPlocs} their locations in the sequence
  whose name is specified by \code{seqname}.

  \code{available.SNPs} returns a character vector containing the names of the
  SNPlocs data packages that are currently available on the Bioconductor
  repositories for your version of R/Bioconductor. A SNPlocs data package
  contains basic SNP information (location and alleles) for a given organism.

  \code{installed.SNPs} returns a character vector containing the names of the
  SNPlocs data packages that are already installed.
}

\note{
  \code{injectSNPs}, \code{SNPlocs_pkgname}, \code{SNPcount} and \code{SNPlocs}
  have the side effect to try to load the SNPlocs data package if it's not
  already loaded.
}

\author{H. Pages}

\seealso{
  \link{BSgenome-class},
  \code{\link[Biostrings]{IUPAC_CODE_MAP}},
  \code{\link[Biostrings]{injectHardMask}},
  \code{\link[Biostrings]{letterFrequencyInSlidingView}},
  \code{\link[Biostrings]{.inplaceReplaceLetterAt}}
}

\examples{
  ## What SNPlocs data packages are already installed:
  installed.SNPs()

  ## What SNPlocs data packages are available:
  available.SNPs()

  if (interactive()) {
    ## Make your choice and install with:
    source("http://bioconductor.org/biocLite.R")
    biocLite("SNPlocs.Hsapiens.dbSNP.20100427")
  }

  ## Inject SNPs from dbSNP into the Human genome:
  library(BSgenome.Hsapiens.UCSC.hg19)
  Hsapiens
  SNPlocs_pkgname(Hsapiens)

  SNP_Hsapiens <- injectSNPs(Hsapiens, "SNPlocs.Hsapiens.dbSNP.20100427")
  SNP_Hsapiens  # note the extra "with SNPs injected from ..." line
  SNPlocs_pkgname(SNP_Hsapiens)
  SNPcount(SNP_Hsapiens)
  head(SNPlocs(SNP_Hsapiens, "chr1"))

  alphabetFrequency(Hsapiens$chr1)
  alphabetFrequency(SNP_Hsapiens$chr1)

  ## Find runs of SNPs of length at least 25 in chr1. Might require
  ## more memory than some platforms can handle (e.g. 32-bit Windows
  ## and maybe some Mac OS X machines with little memory):
  is_32bit_windows <- .Platform$OS.type == "windows" &&
                      .Platform$r_arch == "i386"
  is_macosx <- substr(R.version$os, start=1, stop=6) == "darwin"
  if (!is_32bit_windows && !is_macosx) {
      chr1 <- injectHardMask(SNP_Hsapiens$chr1)
      ambiguous_letters <- paste(DNA_ALPHABET[5:15], collapse="")
      lf <- letterFrequencyInSlidingView(chr1, 25, ambiguous_letters)
      sl <- slice(as.integer(lf), lower=25)
      v1 <- Views(chr1, start(sl), end(sl)+24)
      v1
      max(width(v1))  # length of longest SNP run
  }
}

\keyword{manip}
