\name{XtraSNPlocs-class}
\docType{class}

\alias{class:XtraSNPlocs}
\alias{XtraSNPlocs-class}
\alias{XtraSNPlocs}

% accessors
\alias{provider,XtraSNPlocs-method}
\alias{providerVersion,XtraSNPlocs-method}
\alias{releaseDate,XtraSNPlocs-method}
\alias{releaseName,XtraSNPlocs-method}
\alias{referenceGenome,XtraSNPlocs-method}
\alias{organism,XtraSNPlocs-method}
\alias{species,XtraSNPlocs-method}
\alias{seqinfo,XtraSNPlocs-method}
\alias{seqnames,XtraSNPlocs-method}

% constructor
\alias{newXtraSNPlocs}

% displaying
\alias{show,XtraSNPlocs-method}

% SNP extractors
\alias{snpcount,XtraSNPlocs-method}
\alias{snpsBySeqname}
\alias{snpsBySeqname,XtraSNPlocs-method}
\alias{snpsById}
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

snpsBySeqname(x, seqnames, ...)
\S4method{snpsBySeqname}{XtraSNPlocs}(x, seqnames,
              columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
              drop.rs.prefix=FALSE,
              as.DataFrame=FALSE)

snpsById(x, ids, ...)
\S4method{snpsById}{XtraSNPlocs}(x, ids,
         columns=c("seqnames", "start", "end", "strand", "RefSNP_id"),
         ifnotfound=c("error", "warning", "drop"),
         as.DataFrame=FALSE)
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
    Should the \emph{rs} prefix be dropped from the returned RefSNP ids?
    (RefSNP ids are stored in the \code{RefSNP_id} metadata column of the
    returned object.)
  }
  \item{as.DataFrame}{
    Should the result be returned in a \link[S4Vectors]{DataFrame} instead
    of a \link[GenomicRanges]{GRanges} object?
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
  of the object, except for the columns describing the genomic locations
  of the SNPs (i.e. \code{seqnames}, \code{start}, \code{end}, \code{width},
  and \code{strand}) because they are redundant with the genomic ranges
  represented by the \link[GenomicRanges]{GRanges} object itself.

  Description of the columns:
  \itemize{
    \item \code{seqnames}: The name of the chromosome where each SNP is
          located.
          When a \link[GenomicRanges]{GRanges} object is returned, use
          \code{seqnames()} on it to get this component.

    \item \code{start} and \code{end}: The starting and ending coordinates of
          each SNP with respect to the chromosome indicated in \code{seqnames}.
          Coordinated are 1-based and with respect to the 5' end of the plus
          strand of the chromosome in the reference genome.
          When a \link[GenomicRanges]{GRanges} object is returned, use
          \code{start()}, \code{end()}, or \code{ranges()} on it to get these
          components.

    \item \code{width}: The number of nucleotides spanned by each SNP on the
          reference genome.

    \item \code{strand}: The strand that the alleles of each SNP was reported
          to.

    \item \code{RefSNP_id}: The RefSNP id (a.k.a. \emph{rs id}) of each SNP.

    \item \code{alleles}: The alleles of each SNP.

    \item \code{snpClass}: Class of each SNP. Possible values are
          \code{in-del}, \code{heterozygous}, \code{microsatellite},
          \code{named-locus}, \code{no-variation}, \code{mixed}, and
          \code{multinucleotide-polymorphism}.

    \item \code{loctype}: See \url{ftp://ftp.ncbi.nih.gov/snp/00readme.txt}
          for the 6 loctype codes used by dbSNP and their meanings.
          WARNING: The code assigned to each SNP doesn't seem to be reliable.
          For example, loctypes 1 and 3 officially stand for insertion and
          deletion, but when looking at the data it actually seem to be the
          other way around.
  }
}

\author{H. Pages}

\seealso{
  \itemize{
    \item \code{\link{available.SNPs}}

    \item \link{SNPlocs} objects.
  }
}

\examples{
## COMING SOON!
}

\keyword{methods}
\keyword{classes}
