\name{XtraSNPlocs-class}
\docType{class}

\alias{class:XtraSNPlocs}
\alias{XtraSNPlocs-class}
\alias{XtraSNPlocs}

% getters
\alias{colnames,XtraSNPlocs-method}
\alias{dim,XtraSNPlocs-method}
\alias{provider,XtraSNPlocs-method}
\alias{providerVersion,XtraSNPlocs-method}
\alias{releaseDate,XtraSNPlocs-method}
\alias{releaseName,XtraSNPlocs-method}
\alias{referenceGenome,XtraSNPlocs-method}
\alias{organism,XtraSNPlocs-method}
\alias{commonName,XtraSNPlocs-method}
\alias{species,XtraSNPlocs-method}
\alias{seqinfo,XtraSNPlocs-method}
\alias{seqnames,XtraSNPlocs-method}

% constructor
\alias{newXtraSNPlocs}

% displaying
\alias{show,XtraSNPlocs-method}

% SNP extractors
\alias{snpcount,XtraSNPlocs-method}
\alias{snpsBySeqname,XtraSNPlocs-method}
\alias{snpsByOverlaps,XtraSNPlocs-method}
\alias{snpsById,XtraSNPlocs-method}

\title{XtraSNPlocs objects}

\description{
  The XtraSNPlocs class is a container for storing extra SNP locations and
  alleles for a given organism. 
  While a \link{SNPlocs} object can store only molecular variations of class
  \emph{snp}, an XtraSNPlocs object contains molecular variations of other
  classes (\emph{in-del}, \emph{heterozygous}, \emph{microsatellite},
  \emph{named-locus}, \emph{no-variation}, \emph{mixed},
  \emph{multinucleotide-polymorphism}).

  XtraSNPlocs objects are usually made in advance by a volunteer and made
  available to the Bioconductor community as \emph{XtraSNPlocs data packages}.
  See \code{?\link{available.SNPs}} for how to get the list of
  \emph{SNPlocs and XtraSNPlocs data packages} curently available.

  This man page's main focus is on how to extract data from an
  XtraSNPlocs object.
}

\usage{
\S4method{snpcount}{XtraSNPlocs}(x)

\S4method{snpsBySeqname}{XtraSNPlocs}(x, seqnames,
              columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
              drop.rs.prefix=FALSE,
              as.DataFrame=FALSE)

\S4method{snpsByOverlaps}{XtraSNPlocs}(x, ranges, maxgap=0L, minoverlap=0L,
               type=c("any", "start", "end", "within", "equal"),
               columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
               drop.rs.prefix=FALSE, as.DataFrame=FALSE, ...)

\S4method{snpsById}{XtraSNPlocs}(x, ids,
         columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
         ifnotfound=c("error", "warning", "drop"),
         as.DataFrame=FALSE)

\S4method{colnames}{XtraSNPlocs}(x, do.NULL=TRUE, prefix="col")
}

\arguments{
  \item{x}{
    An XtraSNPlocs object.
  }
  \item{seqnames}{
    The names of the sequences for which to get SNPs. NAs and duplicates
    are not allowed. The supplied \code{seqnames} must be a subset of
    \code{seqlevels(x)}.
  }
  \item{columns}{
    The names of the columns to return. Valid column names are:
    \code{seqnames}, \code{start}, \code{end}, \code{width}, \code{strand},
    \code{RefSNP_id}, \code{alleles}, \code{snpClass}, \code{loctype}.
    See Details section below for a description of these columns.
  }
  \item{drop.rs.prefix}{
    Should the \code{rs} prefix be dropped from the returned RefSNP ids?
    (RefSNP ids are stored in the \code{RefSNP_id} metadata column of the
    returned object.)
  }
  \item{as.DataFrame}{
    Should the result be returned in a \link[S4Vectors]{DataFrame} instead
    of a \link[GenomicRanges]{GRanges} object?
  }
  \item{ranges}{
    One or more regions of interest specified as a \link[GenomicRanges]{GRanges}
    object. A single region of interest can be specified as a character string
    of the form \code{"ch14:5201-5300"}.
  }
  \item{maxgap, minoverlap, type}{
    These arguments are passed to \code{\link[IRanges]{subsetByOverlaps}()}
    which is used internally by \code{snpsByOverlaps}.
    See \code{?IRanges::\link[IRanges]{subsetByOverlaps}} in the \pkg{IRanges}
    package and \code{?GenomicRanges::\link[GenomicRanges]{subsetByOverlaps}} in
    the \pkg{GenomicRanges} package for more information about the
    \code{subsetByOverlaps()} generic and its method for
    \link[GenomicRanges]{GenomicRanges} objects.
  }
  \item{ids}{
    The RefSNP ids to look up (a.k.a. \emph{rs ids}). Can be integer or
    character vector, with or without the \code{"rs"} prefix. NAs are not
    allowed.
  }
  \item{ifnotfound}{
    What to do if SNP ids are not found.
  }
  \item{...}{
    Additional arguments, for use in specific methods.
    Further arguments passed to the \code{snpsByOverlaps} method for
    XtraSNPlocs objects (thru \code{...}) are passed to
    \code{\link[IRanges]{subsetByOverlaps}()}.
  }
  \item{do.NULL, prefix}{
    These arguments are ignored.
  }
}

\value{
  \code{snpcount} returns a named integer vector containing the number
  of SNPs for each chromosome in the reference genome.

  \code{snpsBySeqname} and \code{snpsById} both return a
  \link[GenomicRanges]{GRanges} object with 1 element per SNP,
  unless \code{as.DataFrame} is set to \code{TRUE} in which case they
  return a \link[S4Vectors]{DataFrame} with 1 row per SNP.
  When a \link[GenomicRanges]{GRanges} object is returned, the columns
  requested via the \code{columns} argument are stored as metada columns
  of the object, except for the following columns: \code{seqnames},
  \code{start}, \code{end}, \code{width}, and \code{strand}.
  These "spatial columns" (in the sense that they describe the genomic
  locations of the SNPs) can be accessed by calling the corresponding
  getter on the \link[GenomicRanges]{GRanges} object.

  Summary of available columns (\code{my_snps} being the returned object):
  \itemize{
    \item \code{seqnames}: The name of the chromosome where each SNP is
          located. Access with \code{seqnames(my_snps)} when \code{my_snps}
          is a \link[GenomicRanges]{GRanges} object.

    \item \code{start} and \code{end}: The starting and ending coordinates of
          each SNP with respect to the chromosome indicated in \code{seqnames}.
          Coordinated are 1-based and with respect to the 5' end of the plus
          strand of the chromosome in the reference genome.
          Access with \code{start(my_snps)}, \code{end(my_snps)},
          or \code{ranges(my_snps)} when \code{my_snps} is a
          \link[GenomicRanges]{GRanges} object.

    \item \code{width}: The number of nucleotides spanned by each SNP
          \emph{on the reference genome} (e.g. a width of 0 means the SNP
          is an insertion). Access with \code{width( my_snps)} when
          \code{my_snps} is a \link[GenomicRanges]{GRanges} object.

    \item \code{strand}: The strand that the alleles of each SNP was reported
          to. Access with \code{strand(my_snps)} when \code{my_snps} is a
          \link[GenomicRanges]{GRanges} object.

    \item \code{RefSNP_id}: The RefSNP id (a.k.a. \emph{rs id}) of each SNP.
          Access with \code{mcols(my_snps)$RefSNP_id} when \code{my_snps} is a
          \link[GenomicRanges]{GRanges} object.

    \item \code{alleles}: The alleles of each SNP in the format used by dbSNP.
          Access with \code{mcols(my_snps)$alleles} when \code{my_snps} is a
          \link[GenomicRanges]{GRanges} object.

    \item \code{snpClass}: Class of each SNP. Possible values are
          \code{in-del}, \code{heterozygous}, \code{microsatellite},
          \code{named-locus}, \code{no-variation}, \code{mixed}, and
          \code{multinucleotide-polymorphism}.
          Access with \code{mcols(my_snps)$snpClass} when \code{my_snps} is a
          \link[GenomicRanges]{GRanges} object.

    \item \code{loctype}: See \url{ftp://ftp.ncbi.nih.gov/snp/00readme.txt}
          for the 6 loctype codes used by dbSNP, and their meanings.
          WARNING: The code assigned to each SNP doesn't seem to be reliable.
          For example, loctype codes 1 and 3 officially stand for insertion
          and deletion, respectively. However, when looking at the SNP ranges
          it actually seems to be the other way around.
          Access with \code{mcols(my_snps)$loctype} when \code{my_snps} is a
          \link[GenomicRanges]{GRanges} object.
  }

  \code{colnames(x)} returns the names of the available columns.
}

\author{H. Pages}

\seealso{
  \itemize{
    \item \code{\link{available.SNPs}}

    \item \link{SNPlocs} objects.
  }
}

\examples{
library(XtraSNPlocs.Hsapiens.dbSNP141.GRCh38)
snps <- XtraSNPlocs.Hsapiens.dbSNP141.GRCh38
snpcount(snps)
colnames(snps)

## ---------------------------------------------------------------------
## snpsBySeqname()
## ---------------------------------------------------------------------
## Get the location, RefSNP id, and alleles for all "extra SNPs"
## located on chromosome 22 and MT:
snpsBySeqname(snps, c("ch22", "chMT"), columns=c("RefSNP_id", "alleles"))

## ---------------------------------------------------------------------
## snpsByOverlaps()
## ---------------------------------------------------------------------
## Get the location, RefSNP id, and alleles for all "extra SNPs"
## overlapping some regions of interest:
snpsByOverlaps(snps, "ch22:33.63e6-33.64e6",
               columns=c("RefSNP_id", "alleles"))

## With the regions of interest being all the known CDS for hg38
## (except for the chromosome naming convention, hg38 is the same
## as GRCh38):
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
hg38_cds <- cds(txdb)
seqlevelsStyle(hg38_cds)  # UCSC
seqlevelsStyle(snps)  # dbSNP
seqlevelsStyle(hg38_cds) <- seqlevelsStyle(snps)
genome(hg38_cds) <- genome(snps)
snpsByOverlaps(snps, hg38_cds, columns=c("RefSNP_id", "alleles"))

## ---------------------------------------------------------------------
## snpsById()
## ---------------------------------------------------------------------
## Get the location and alleles for some RefSNP ids:
my_rsids <- c("rs367617508", "rs398104919", "rs3831697", "rs372470289",
              "rs141568169", "rs34628976", "rs67551854")
snpsById(snps, my_rsids, c("RefSNP_id", "alleles"))

## See ?XtraSNPlocs.Hsapiens.dbSNP141.GRCh38 for examples of using
## snpsBySeqname() and snpsById().
}

\keyword{methods}
\keyword{classes}
