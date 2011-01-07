\name{BSgenomeForge}

\alias{BSgenomeForge}

\alias{BSgenomeDataPkgSeed}
\alias{class:BSgenomeDataPkgSeed}
\alias{BSgenomeDataPkgSeed-class}

\alias{forgeBSgenomeDataPkg}
\alias{forgeBSgenomeDataPkg,BSgenomeDataPkgSeed-method}
\alias{forgeBSgenomeDataPkg,list-method}
\alias{forgeBSgenomeDataPkg,character-method}
\alias{forgeSeqlengthsFile}
\alias{forgeSeqFiles}
\alias{forgeMasksFiles}


\title{The BSgenomeForge functions}

\description{
  A set of functions for making a BSgenome data package.
}

\usage{
## Top-level BSgenomeForge function:

forgeBSgenomeDataPkg(x, seqs_srcdir=".", masks_srcdir=".", destdir=".", verbose=TRUE)

## Low-level BSgenomeForge functions:

forgeSeqlengthsFile(seqnames, prefix="", suffix=".fa",
                    seqs_srcdir=".", seqs_destdir=".", verbose=TRUE)

forgeSeqFiles(seqnames, mseqnames=NULL, prefix="", suffix=".fa",
              seqs_srcdir=".", seqs_destdir=".", verbose=TRUE)

forgeMasksFiles(seqnames, nmask_per_seq,
                seqs_destdir=".", masks_srcdir=".", masks_destdir=".",
                AGAPSfiles_type="gap", AGAPSfiles_name=NA,
                AGAPSfiles_prefix="", AGAPSfiles_suffix="_gap.txt",
                RMfiles_name=NA, RMfiles_prefix="", RMfiles_suffix=".fa.out",
                TRFfiles_name=NA, TRFfiles_prefix="", TRFfiles_suffix=".bed",
                verbose=TRUE)
}

\arguments{
  \item{x}{
    A BSgenomeDataPkgSeed object or the name of a BSgenome data package seed
    file.
    See the BSgenomeForge vignette in this package for more information.
  }
  \item{seqs_srcdir, masks_srcdir}{
    Single strings indicating the path to the source directories i.e. to the
    directories containing the source data files.
    Only read access to these directories is needed.
    See the BSgenomeForge vignette in this package for more information.
  }
  \item{destdir}{
    A single string indicating the path to the directory where the source tree
    of the target package should be created.
    This directory must already exist.
    See the BSgenomeForge vignette in this package for more information.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{seqnames, mseqnames}{
    A character vector containing the names of the single (for \code{seqnames})
    and multiple (for \code{mseqnames}) sequences to forge.
    See the BSgenomeForge vignette in this package for more information.
  }
  \item{prefix, suffix}{
    See the BSgenomeForge vignette in this package for more information,
    in particular the description of the \code{seqfiles_prefix} and
    \code{seqfiles_suffix} fields of a BSgenome data package seed file.
  }
  \item{seqs_destdir, masks_destdir}{
    During the forging process the source data files are converted into
    serialized Biostrings objects. \code{seqs_destdir} and \code{masks_destdir}
    must be single strings indicating the path to the directories where these
    serialized objects should be saved.
    These directories must already exist.

    \code{forgeSeqlengthsFile} will produce a single .rda file.
    Both \code{forgeSeqFiles} and \code{forgeMasksFiles} will produce one
    .rda file per sequence.
  }
  \item{nmask_per_seq}{
    A single integer indicating the desired number of masks per sequence.
    See the BSgenomeForge vignette in this package for more information.
  }
  \item{AGAPSfiles_type, AGAPSfiles_name, AGAPSfiles_prefix, AGAPSfiles_suffix,
        RMfiles_name, RMfiles_prefix, RMfiles_suffix,
        TRFfiles_name, TRFfiles_prefix, TRFfiles_suffix}{
    These arguments are named accordingly to the corresponding fields of a
    BSgenome data package seed file.
    See the BSgenomeForge vignette in this package for more information.
  }
}

\details{
  These functions are intended for Bioconductor users who want to make a new
  BSgenome data package, not for regular users of these packages.
  See the BSgenomeForge vignette in this package
  (\code{vignette("BSgenomeForge")}) for an extensive coverage
  of this topic.
}

\author{H. Pages}

\examples{
  forgeSeqFiles("chrM", prefix="ce2", suffix=".fa",
                seqs_srcdir=system.file("extdata", package="BSgenome"),
                seqs_destdir=tempdir())
  load(file.path(tempdir(), "chrM.rda"))
  chrM
}

\keyword{manip}
