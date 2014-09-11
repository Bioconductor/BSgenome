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
  \item{...}{
    Additional arguments, for use in specific methods.
  }
}

\value{
  \code{snpcount} returns a named integer vector containing the number
  of SNPs for each sequence in the reference genome.

  By default (i.e. when \code{as.GRanges=FALSE}), \code{snplocs} returns a
  data frame with 1 row per SNP and the following columns:
  \enumerate{
    \item \code{RefSNP_id}: RefSNP ID (aka "rs id") with \code{"rs"}
          prefix removed. Character vector with no NAs and no duplicates.
    \item \code{alleles_as_ambig}: A character vector with no NAs
          containing the alleles for each SNP represented by an IUPAC
          nucleotide ambiguity code.
          See \code{?\link[Biostrings]{IUPAC_CODE_MAP}} in the
          \pkg{Biostrings} package for more information.
    \item \code{loc}: The 1-based location of the SNP relative to the
          first base at the 5' end of the plus strand of the reference
          sequence.
  }
  Otherwise (i.e. when \code{as.GRanges=TRUE}), it returns a
  \link[GenomicRanges]{GRanges} object with extra columns
  \code{"RefSNP_id"} and \code{"alleles_as_ambig"}.
  Note that all the elements (genomic ranges) in this
  \link[GenomicRanges]{GRanges} object have their strand set
  to \code{"+"} and that all the sequence lengths are set to \code{NA}.

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
## COMING SOON!
}

\keyword{methods}
\keyword{classes}
