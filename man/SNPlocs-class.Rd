\name{SNPlocs-class}
\docType{class}

\alias{class:SNPlocs}
\alias{SNPlocs-class}
\alias{SNPlocs}

\alias{class:ODLT_SNPlocs}
\alias{ODLT_SNPlocs-class}
\alias{ODLT_SNPlocs}

\alias{class:OldFashionSNPlocs}
\alias{OldFashionSNPlocs-class}
\alias{OldFashionSNPlocs}

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
\alias{seqinfo,SNPlocs-method}
\alias{seqnames,SNPlocs-method}

% constructors
\alias{new_ODLT_SNPlocs}
\alias{newSNPlocs}

% displaying
\alias{show,SNPlocs-method}

% SNP extractors
\alias{snpcount}
\alias{snpcount,SNPlocs-method}
\alias{snpcount,ODLT_SNPlocs-method}
\alias{snpcount,OldFashionSNPlocs-method}
\alias{snplocs}
\alias{snplocs,SNPlocs-method}
\alias{snplocs,ODLT_SNPlocs-method}
\alias{snplocs,OldFashionSNPlocs-method}
\alias{snpsBySeqname}
\alias{snpsBySeqname,SNPlocs-method}
\alias{snpsBySeqname,ODLT_SNPlocs-method}
\alias{snpsBySeqname,OldFashionSNPlocs-method}
\alias{snpsByOverlaps}
\alias{snpsByOverlaps,SNPlocs-method}
\alias{snpsByOverlaps,ODLT_SNPlocs-method}
\alias{snpsByOverlaps,OldFashionSNPlocs-method}
\alias{snpsById}
\alias{snpsById,SNPlocs-method}
\alias{snpsById,ODLT_SNPlocs-method}
\alias{snpsById,OldFashionSNPlocs-method}

% Old SNPlocs extractors (deprecated)
\alias{snpid2loc}
\alias{snpid2loc,SNPlocs-method}
\alias{snpid2loc,OldFashionSNPlocs-method}
\alias{snpid2alleles}
\alias{snpid2alleles,SNPlocs-method}
\alias{snpid2alleles,OldFashionSNPlocs-method}
\alias{snpid2grange}
\alias{snpid2grange,SNPlocs-method}
\alias{snpid2grange,OldFashionSNPlocs-method}

\title{SNPlocs objects}

\description{
  The SNPlocs class is a container for storing known SNP locations (of class
  \emph{snp}) for a given organism.

  SNPlocs objects are usually made in advance by a volunteer and made
  available to the Bioconductor community as \emph{SNPlocs data packages}.
  See \code{?\link{available.SNPs}} for how to get the list of
  \emph{SNPlocs and \link{XtraSNPlocs} data packages} curently available.

  The main focus of this man page is on how to extract SNPs from an
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
}

\arguments{
  \item{x}{
    A SNPlocs object.
  }
  \item{seqnames}{
    The names of the sequences for which to get SNPs. Must be a subset of
    \code{seqlevels(x)}. NAs and duplicates are not allowed.
  }
  \item{...}{
    Additional arguments, for use in specific methods.

    Arguments passed to the \code{snpsByOverlaps} method for SNPlocs
    objects thru \code{...} are passed to internal call to
    \code{\link[IRanges]{subsetByOverlaps}()}.
  }
  \item{drop.rs.prefix}{
    Should the \code{rs} prefix be dropped from the returned RefSNP ids?
    (RefSNP ids are stored in the \code{RefSNP_id} metadata column of the
    returned object.)
  }
  \item{ranges}{
    One or more genomic regions of interest specified as a
    \link[GenomicRanges]{GRanges} or \link[GenomicRanges]{GPos} object.
    A single region of interest can be specified as a character string of
    the form \code{"ch14:5201-5300"}.
  }
  \item{maxgap, minoverlap, type}{
    These arguments are passed to \code{\link[IRanges]{subsetByOverlaps}()}
    which is used internally by \code{snpsByOverlaps}.

    Note that the default value for \code{minoverlap} is 0 which means
    that, by default, in addition to the SNPs that are located \emph{within}
    the genomic regions specified thru the \code{ranges} argument,
    \code{snpsByOverlaps} also returns SNPs that are \emph{adjacent} to
    these regions. Use \code{minoverlap=1L} to omit these SNPs.
  }
  \item{ids}{
    The RefSNP ids to look up (a.k.a. rs ids). Can be integer or character
    vector, with or without the \code{"rs"} prefix. NAs are not allowed.
  }
  \item{ifnotfound}{
    What to do if SNP ids are not found.
  }
}

\value{
  \code{snpcount} returns a named integer vector containing the number
  of SNPs for each sequence in the reference genome.

  \code{snpsBySeqname}, \code{snpsByOverlaps}, and \code{snpsById} return
  an \emph{unstranded} \link[GenomicRanges]{GPos} object with 1 element
  (genomic position) per SNP and the following metadata columns:
  \itemize{
    \item \code{RefSNP_id}: RefSNP ID (aka "rs id"). Character vector
          with no NAs and no duplicates.
    \item \code{alleles_as_ambig}: A character vector with no NAs
          containing the alleles for each SNP represented by an IUPAC
          nucleotide ambiguity code.
          See \code{?\link[Biostrings]{IUPAC_CODE_MAP}} in the
          \pkg{Biostrings} package for more information.
          The alleles are always reported with respect to the positive
          strand.
  }
  Note that this \link[GenomicRanges]{GPos} object is \emph{unstranded} i.e.
  all the SNPs in it have their strand set to \code{"*"}.

  If \code{ifnotfound="error"}, the object returned by \code{snpsById}
  is guaranteed to be \emph{parallel} to \code{ids}, that is, the i-th
  element in the \link[GenomicRanges]{GPos} object corresponds to the
  i-th element in \code{ids}.
}

\author{H. Pagès}

\seealso{
  \itemize{
    \item \code{\link{available.SNPs}}

    \item \link[GenomicRanges]{GPos} and \link[GenomicRanges]{GRanges}
          objects in the \pkg{GenomicRanges} package.

    \item \link{XtraSNPlocs} packages and objects for molecular variations
          of class other than \emph{snp} e.g. of class \emph{in-del},
          \emph{heterozygous}, \emph{microsatellite}, etc...

    \item \code{IRanges::\link[IRanges]{subsetByOverlaps}} in the
          \pkg{IRanges} package and
          \code{GenomicRanges::\link[GenomicRanges]{subsetByOverlaps}}
          in the \pkg{GenomicRanges} package for more information about the
          \code{subsetByOverlaps()} generic and its method for
          \link[GenomicRanges]{GenomicRanges} objects.

    \item \code{\link{injectSNPs}}

    \item \code{\link[Biostrings]{IUPAC_CODE_MAP}} in the \pkg{Biostrings}
          package.
  }
}

\examples{
library(SNPlocs.Hsapiens.dbSNP144.GRCh38)
snps <- SNPlocs.Hsapiens.dbSNP144.GRCh38
snpcount(snps)

## ---------------------------------------------------------------------
## snpsBySeqname()
## ---------------------------------------------------------------------

## Get all SNPs located on chromosome 22 or MT:
snpsBySeqname(snps, c("22", "MT"))

## ---------------------------------------------------------------------
## snpsByOverlaps()
## ---------------------------------------------------------------------

## Get all SNPs overlapping some genomic region of interest:
snpsByOverlaps(snps, "22:33.63e6-33.64e6")

## With the regions of interest being all the known CDS for hg38
## located on chromosome 22 or MT (except for the chromosome naming
## convention, hg38 is the same as GRCh38):
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
my_cds <- cds(txdb)
seqlevels(my_cds, pruning.mode="coarse") <- c("chr22", "chrM")
seqlevelsStyle(my_cds)  # UCSC
seqlevelsStyle(snps)    # NCBI
seqlevelsStyle(my_cds) <- seqlevelsStyle(snps)
genome(my_cds) <- genome(snps)
my_snps <- snpsByOverlaps(snps, my_cds)
my_snps
table(my_snps \%within\% my_cds)

## Note that, by default, in addition to the SNPs that are located
## within the regions of interest, snpsByOverlaps() also returns SNPs
## that are adjacent to these regions. Use 'minoverlap=1L' to omit these
## SNPs:
my_snps2 <- snpsByOverlaps(snps, my_cds, minoverlap=1)
table(my_snps2 \%within\% my_cds)

## SNPs adjacent to the regions of interest:
setdiff(my_cds, my_snps2)

## ---------------------------------------------------------------------
## snpsById()
## ---------------------------------------------------------------------

## Lookup some RefSNP ids:
my_rsids <- c("rs10458597", "rs12565286", "rs7553394")
\dontrun{
  snpsById(snps, my_rsids)  # error, rs7553394 not found
}
## This example requires too much memory on 32-bit Windows!
if (.Platform$OS.type != "windows" || .Platform$r_arch != "i386") {
    snpsById(snps, my_rsids, ifnotfound="drop")
}
}

\keyword{methods}
\keyword{classes}
