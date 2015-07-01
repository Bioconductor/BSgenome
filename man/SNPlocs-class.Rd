\name{SNPlocs-class}
\docType{class}

\alias{class:SNPlocs}
\alias{SNPlocs-class}
\alias{SNPlocs}

% accessors
\alias{provider,SNPlocs-method}
\alias{providerVersion,SNPlocs-method}
\alias{releaseDate,SNPlocs-method}
\alias{releaseName,SNPlocs-method}
\alias{referenceGenome}
\alias{referenceGenome,SNPlocs-method}
\alias{compatibleGenomes}
\alias{compatibleGenomes,SNPlocs-method}
\alias{organism,SNPlocs-method}
\alias{commonName,SNPlocs-method}
\alias{species,SNPlocs-method}
\alias{seqinfo,SNPlocs-method}
\alias{seqnames,SNPlocs-method}

% constructor
\alias{newSNPlocs}

% displaying
\alias{show,SNPlocs-method}

% SNP extractors
\alias{snpcount}
\alias{snpcount,SNPlocs-method}
\alias{snpsBySeqname}
\alias{snpsBySeqname,SNPlocs-method}
\alias{snpsByOverlaps}
\alias{snpsByOverlaps,SNPlocs-method}
\alias{snpsById}
\alias{snpsById,SNPlocs-method}

% OLD API
\alias{snplocs}
\alias{snplocs,SNPlocs-method}
\alias{snpid2loc}
\alias{snpid2loc,SNPlocs-method}
\alias{snpid2alleles}
\alias{snpid2alleles,SNPlocs-method}
\alias{snpid2grange}
\alias{snpid2grange,SNPlocs-method}

\title{SNPlocs objects}

\description{
  The SNPlocs class is a container for storing known SNP locations for a
  given organism. SNPlocs objects are usually made in advance by
  a volunteer and made available to the Bioconductor community as
  "SNPlocs data packages".
  See \code{?\link{available.SNPs}} for how to get the list of
  "SNPlocs data packages" curently available.

  This man page's main focus is on how to extract information from a
  SNPlocs object.
}

\usage{
snpcount(x)

snpsBySeqname(x, seqnames, ...)

\S4method{snpsBySeqname}{SNPlocs}(x, seqnames, drop.rs.prefix=FALSE)

snpsByOverlaps(x, ranges, maxgap=0L, minoverlap=0L,
               type=c("any", "start", "end", "within", "equal"), ...)

\S4method{snpsByOverlaps}{SNPlocs}(x, ranges, maxgap=0L, minoverlap=0L,
               type=c("any", "start", "end", "within", "equal"),
               drop.rs.prefix=FALSE, ...)

snpsById(x, ids, ...)

\S4method{snpsById}{SNPlocs}(x, ids, ifnotfound=c("error", "warning", "drop"))

## Old API
## ------------------------------------

snplocs(x, seqname, ...)
\S4method{snplocs}{SNPlocs}(x, seqname, as.GRanges=FALSE, caching=TRUE)

snpid2loc(x, snpid, ...)
\S4method{snpid2loc}{SNPlocs}(x, snpid, caching=TRUE)

snpid2alleles(x, snpid, ...)
\S4method{snpid2alleles}{SNPlocs}(x, snpid, caching=TRUE)

snpid2grange(x, snpid, ...)
\S4method{snpid2grange}{SNPlocs}(x, snpid, caching=TRUE)
}

\arguments{
  \item{x}{
    A SNPlocs object.
  }
  \item{seqnames}{
    The names of the sequences for which to get SNPs. NAs and duplicates
    are not allowed. The supplied \code{seqnames} must be a subset of
    \code{seqlevels(x)}.
  }
  \item{...}{
    Additional arguments, for use in specific methods.

    Further arguments passed to the \code{snpsByOverlaps} method
    for SNPlocs objects (thru \code{...}) are passed to
    \code{\link[IRanges]{subsetByOverlaps}()}.
  }
  \item{drop.rs.prefix}{
    Should the \code{rs} prefix be dropped from the returned RefSNP ids?
    (RefSNP ids are stored in the \code{RefSNP_id} metadata column of the
    returned object.)
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
  \item{seqname}{
    The name of the sequence for which to get the SNP locations
    and alleles.

    If \code{as.GRanges} is \code{FALSE}, only one sequence can
    be specified (i.e. \code{seqname} must be a single string).
    If \code{as.GRanges} is \code{TRUE}, an arbitrary number of
    sequences can be specified (i.e. \code{seqname} can be
    a character vector of arbitrary length).
  }
  \item{as.GRanges}{
    \code{TRUE} or \code{FALSE}. If \code{TRUE}, then the SNP locations
    and alleles are returned in a \link[GenomicRanges]{GRanges} object.
    Otherwise (the default), they are returned in a data frame (see below).
  }
  \item{caching}{
    Should the loaded SNPs be cached in memory for faster further
    retrieval but at the cost of increased memory usage?
  }
  \item{snpid}{
    The SNP ids to look up (e.g. rs ids). Can be integer or character vector,
    with or without the \code{"rs"} prefix. NAs are not allowed.
  }
}

\value{
  \code{snpcount} returns a named integer vector containing the number
  of SNPs for each sequence in the reference genome.

  \code{snpsBySeqname}, \code{snpsByOverlaps}, and \code{snpsById} return
  a \link[GenomicRanges]{GRanges} object with 1 element (genomic range)
  per SNP and the following metadata columns:
  \itemize{
    \item \code{RefSNP_id}: RefSNP ID (aka "rs id"). Character vector
          with no NAs and no duplicates.
    \item \code{alleles_as_ambig}: A character vector with no NAs
          containing the alleles for each SNP represented by an IUPAC
          nucleotide ambiguity code.
          See \code{?\link[Biostrings]{IUPAC_CODE_MAP}} in the
          \pkg{Biostrings} package for more information.
  }
  Note that all the elements (genomic ranges) in this
  \link[GenomicRanges]{GRanges} object have their strand set to \code{"+"}.

  If \code{ifnotfound="error"}, the object returned by \code{snpsById}
  is guaranteed to be \emph{parallel} to \code{ids}, that is, the i-th
  element in the \link[GenomicRanges]{GRanges} object corresponds to the
  i-th element in \code{ids}.

  \subsection{Old API}{
    Note that \code{snplocs} is superseded by \code{snpsBySeqname}, and
    \code{snpid2loc}, \code{snpid2alleles}, and \code{snpid2grange} are
    superseded by \code{snpsById}.

    By default (i.e. when \code{as.GRanges=FALSE}), \code{snplocs} returns a
    data frame with 1 row per SNP and the following columns:
    \enumerate{
      \item \code{RefSNP_id}: Same as above but with \code{"rs"} prefix
            always removed.
      \item \code{alleles_as_ambig}: Same as above.
      \item \code{loc}: The 1-based location of the SNP relative to the
            first base at the 5' end of the plus strand of the reference
            sequence.
    }
    Otherwise (i.e. when \code{as.GRanges=TRUE}), it returns a
    \link[GenomicRanges]{GRanges} object with metadata columns
    \code{"RefSNP_id"} and \code{"alleles_as_ambig"}.

    \code{snpid2loc} and \code{snpid2alleles} both return a named vector
    (integer vector for the former, character vector for the latter)
    where each (name, value) pair corresponds to a supplied SNP id.
    For both functions the name in (name, value) is the chromosome
    of the SNP id. The value in (name, value) is the position of the
    SNP id on the chromosome for \code{snpid2loc}, and a single IUPAC
    code representing the associated alleles for \code{snpid2alleles}.

    \code{snpid2grange} returns a \link[GenomicRanges]{GRanges} object
    similar to the one returned by \code{snplocs} (when used with
    \code{as.GRanges=TRUE}) and where each element corresponds to a
    supplied SNP id.
  }
}

\author{H. Pages}

\seealso{
  \itemize{
    \item \code{\link{available.SNPs}}

    \item \code{\link{injectSNPs}}

    \item \code{\link[Biostrings]{IUPAC_CODE_MAP}} in the \pkg{Biostrings}
          package.
  }
}

\examples{
library(SNPlocs.Hsapiens.dbSNP141.GRCh38)
snps <- SNPlocs.Hsapiens.dbSNP141.GRCh38
snpcount(snps)

## ---------------------------------------------------------------------
## snpsBySeqname()
## ---------------------------------------------------------------------
## Get all SNPs located on chromosome 22 and MT:
snpsBySeqname(snps, c("ch22", "chMT"))

## ---------------------------------------------------------------------
## snpsByOverlaps()
## ---------------------------------------------------------------------
## Get all SNPs overlapping some regions of interest:
snpsByOverlaps(snps, "ch22:33.63e6-33.64e6")

## With the regions of interest being all the known CDS for hg38
## located on chr22 or chrMT (except for the chromosome naming
## convention, hg38 is the same as GRCh38):
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
my_cds <- cds(txdb)
seqlevels(my_cds, force=TRUE) <- c("chr22", "chrMT")
seqlevelsStyle(my_cds)  # UCSC
seqlevelsStyle(snps)  # dbSNP
seqlevelsStyle(my_cds) <- seqlevelsStyle(snps)
genome(my_cds) <- genome(snps)
snpsByOverlaps(snps, my_cds)

## ---------------------------------------------------------------------
## snpsById()
## ---------------------------------------------------------------------
## Lookup some RefSNP ids:
my_rsids <- c("rs10458597", "rs12565286", "rs7553394")
\dontrun{
  snpsById(snps, my_rsids)  # error, rs7553394 not found
}
snpsById(snps, my_rsids, ifnotfound="drop")
}

\keyword{methods}
\keyword{classes}
