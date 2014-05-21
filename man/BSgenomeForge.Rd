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

forgeBSgenomeDataPkg(x, seqs_srcdir=".", destdir=".", verbose=TRUE)

## Low-level BSgenomeForge functions:

forgeSeqlengthsFile(seqnames, prefix="", suffix=".fa",
                    seqs_srcdir=".", seqs_destdir=".", verbose=TRUE)

forgeSeqFiles(seqnames, mseqnames=NULL,
              seqfile_name=NA, prefix="", suffix=".fa",
              seqs_srcdir=".", seqs_destdir=".",
              ondisk_seq_format=c("2bit", "rda", "fa.rz", "fa"),
              verbose=TRUE)

forgeMasksFiles(seqnames, nmask_per_seq,
                seqs_destdir=".",
                ondisk_seq_format=c("2bit", "rda", "fa.rz", "fa"),
                masks_srcdir=".", masks_destdir=".",
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
  \item{ondisk_seq_format}{
    Specifies how the single sequences should be stored in the forged package.
    Can be \code{"2bit"}, \code{"rda"}, \code{"fa.rz"}, or \code{"fa"}.
    If \code{"2bit"} (the default), then all the single sequences are stored
    in a single twoBit file.
    If \code{"rda"}, then each single sequence is stored in a separated
    serialized \link[Biostrings]{XString} object (one per single sequence).
    If \code{"fa.rz"} or \code{"fa"}, then all the single sequences are stored
    in a single FASTA file (compressed in the RAZip format if \code{"fa.rz"}).
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{seqnames, mseqnames}{
    A character vector containing the names of the single (for \code{seqnames})
    and multiple (for \code{mseqnames}) sequences to forge.
    See the BSgenomeForge vignette in this package for more information.
  }
  \item{seqfile_name, prefix, suffix}{
    See the BSgenomeForge vignette in this package for more information,
    in particular the description of the \code{seqfile_name},
    \code{seqfiles_prefix} and \code{seqfiles_suffix} fields of a BSgenome
    data package seed file.
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
seqs_srcdir <- system.file("extdata", package="BSgenome")
seqnames <- c("chrX", "chrM")

## Forge .rda sequence files:
forgeSeqFiles(seqnames, prefix="ce2", suffix=".fa.gz",
              seqs_srcdir=seqs_srcdir,
              seqs_destdir=tempdir(), ondisk_seq_format="rda")

## Forge .2bit sequence files:
forgeSeqFiles(seqnames, prefix="ce2", suffix=".fa.gz",
              seqs_srcdir=seqs_srcdir,
              seqs_destdir=tempdir(), ondisk_seq_format="2bit")

## Sanity checks:
library(BSgenome.Celegans.UCSC.ce2)
genome <- BSgenome.Celegans.UCSC.ce2

load(file.path(tempdir(), "chrX.rda"))
stopifnot(genome$chrX == chrX)
load(file.path(tempdir(), "chrM.rda"))
stopifnot(genome$chrM == chrM)

ce2_sequences <- import(file.path(tempdir(), "single_sequences.2bit"))
ce2_sequences0 <- DNAStringSet(list(chrX=genome$chrX, chrM=genome$chrM))
stopifnot(identical(names(ce2_sequences0), names(ce2_sequences)) &&
          all(ce2_sequences0 == ce2_sequences))
}

\keyword{manip}
