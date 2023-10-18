\name{BSgenomeForge}

\alias{BSgenomeForge}

\alias{BSgenomeDataPkgSeed}
\alias{class:BSgenomeDataPkgSeed}
\alias{BSgenomeDataPkgSeed-class}

\alias{forgeBSgenomeDataPkg}
\alias{forgeBSgenomeDataPkg,BSgenomeDataPkgSeed-method}
\alias{forgeBSgenomeDataPkg,list-method}
\alias{forgeBSgenomeDataPkg,character-method}
\alias{forgeSeqlengthsRdsFile}
\alias{forgeSeqlengthsRdaFile}
\alias{forgeSeqFiles}
\alias{forgeMasksFiles}


\title{The BSgenomeForge functions}

\description{
  A set of functions for making a BSgenome data package.

  IMPORTANT NOTE: A new package, the \pkg{BSgenomeForge} package,
  provides more user-friendly tools for creating a BSgenome data package
  from an NCBI assembly or UCSC genome. However, if your assembly or genome
  is not from NCBI or UCSC, you should still use the tools documented below.
}

\usage{
## Top-level BSgenomeForge function:

forgeBSgenomeDataPkg(x, seqs_srcdir=".", destdir=".", replace=FALSE, verbose=TRUE)

## Low-level BSgenomeForge functions:

forgeSeqlengthsRdsFile(seqnames, prefix="", suffix=".fa",
                       seqs_srcdir=".", seqs_destdir=".",
                       genome=NA_character_, verbose=TRUE)

forgeSeqlengthsRdaFile(seqnames, prefix="", suffix=".fa",
                       seqs_srcdir=".", seqs_destdir=".",
                       genome=NA_character_, verbose=TRUE)

forgeSeqFiles(provider, genome,
              seqnames, mseqnames=NULL,
              seqfile_name=NA, prefix="", suffix=".fa",
              seqs_srcdir=".", seqs_destdir=".",
              ondisk_seq_format=c("2bit", "rds", "rda", "fa.rz", "fa"),
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
  \item{replace}{
    \code{TRUE} or \code{FALSE}.
    When set to TRUE, \code{replace} replaces the package directory
    if it already exists.
  }
  \item{verbose}{
    \code{TRUE} or \code{FALSE}.
  }
  \item{provider}{
    The provider of the \emph{sequence data files} e.g.
    \code{"UCSC"}, \code{"NCBI"}, \code{"BDGP"}, \code{"FlyBase"}, etc...
  }
  \item{genome}{
    The name of the genome. Typically the name of an NCBI assembly (e.g.
    \code{"GRCh38.p12"}, \code{"WBcel235"}, \code{"TAIR10.1"},
    \code{"ARS-UCD1.2"}, etc...) or UCSC genome (e.g. \code{"hg38"},
    \code{"bosTau9"}, \code{"galGal6"}, \code{"ce11"}, etc...).
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

    Both \code{forgeSeqlengthsRdsFile} and \code{forgeSeqlengthsRdaFile}
    will produce a single \code{.rds} or \code{.rda} file.
    Both \code{forgeSeqFiles} and \code{forgeMasksFiles} will produce one
    file per sequence (all files being either \code{.rds} or \code{.rda}
    files).
  }
  \item{ondisk_seq_format}{
    Specifies how the single sequences should be stored in the forged package.
    Can be \code{"2bit"}, \code{"rds"}, \code{"rda"}, \code{"fa.rz"},
    or \code{"fa"}.
    If \code{"2bit"} (the default), then all the single sequences are stored
    in a single twoBit file.
    If \code{"rds"} or \code{"rda"}, then each single sequence is stored in
    a separated serialized \link[Biostrings]{XString} derivative (one per
    single sequence).
    If \code{"fa.rz"} or \code{"fa"}, then all the single sequences are stored
    in a single FASTA file (compressed in the RAZip format if \code{"fa.rz"}).
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

\author{H. Pagès}

\seealso{
  \itemize{
    \item \code{\link[BSgenomeForge]{forgeBSgenomeDataPkgFromNCBI}} and
          \code{\link[BSgenomeForge]{forgeBSgenomeDataPkgFromUCSC}} in the
          \pkg{BSgenomeForge} package.

    \item \code{\link{available.genomes}} to find BSgenome data packages
          available in Bioconductor.

    \item \link{BSgenome} objects.
  }
}

\examples{
seqs_srcdir <- system.file("extdata", package="BSgenome")
seqnames <- c("chrX", "chrM")

## Forge .2bit sequence files:
forgeSeqFiles("UCSC", "ce2",
              seqnames, prefix="ce2", suffix=".fa.gz",
              seqs_srcdir=seqs_srcdir,
              seqs_destdir=tempdir(), ondisk_seq_format="2bit")

## Forge .rds sequence files:
forgeSeqFiles("UCSC", "ce2",
              seqnames, prefix="ce2", suffix=".fa.gz",
              seqs_srcdir=seqs_srcdir,
              seqs_destdir=tempdir(), ondisk_seq_format="rds")

## Sanity checks:
library(BSgenome.Celegans.UCSC.ce2)
genome <- BSgenome.Celegans.UCSC.ce2

ce2_sequences <- import(file.path(tempdir(), "single_sequences.2bit"))
ce2_sequences0 <- DNAStringSet(list(chrX=genome$chrX, chrM=genome$chrM))
stopifnot(identical(names(ce2_sequences0), names(ce2_sequences)),
          all(ce2_sequences0 == ce2_sequences))

chrX <- readRDS(file.path(tempdir(), "chrX.rds"))
stopifnot(genome$chrX == chrX)
chrM <- readRDS(file.path(tempdir(), "chrM.rds"))
stopifnot(genome$chrM == chrM)
}

\keyword{manip}
