\name{GenomicRanges-class}
\docType{class}

% Class
\alias{class:GenomicRanges}
\alias{GenomicRanges-class}
\alias{GenomicRanges}

% Constructor-like function:
\alias{GenomicRanges}

% Coercion methods:
\alias{coerce,RangedData,GenomicRanges-method}
\alias{as.data.frame,GenomicRanges-method}

% Accessor methods:
\alias{seqnames,GenomicRanges-method}
\alias{seqnames<-,GenomicRanges-method}
\alias{ranges,GenomicRanges-method}
\alias{ranges<-,GenomicRanges-method}
\alias{strand,GenomicRanges-method}
\alias{strand<-,GenomicRanges-method}
\alias{values,GenomicRanges-method}
\alias{values<-,GenomicRanges-method}
\alias{names,GenomicRanges-method}
\alias{names<-,GenomicRanges-method}

% Ranges methods:
\alias{start,GenomicRanges-method}
\alias{start<-,GenomicRanges-method}
\alias{end,GenomicRanges-method}
\alias{end<-,GenomicRanges-method}
\alias{width,GenomicRanges-method}
\alias{width<-,GenomicRanges-method}

% DataTable methods:
\alias{[,GenomicRanges-method}
\alias{[<-,GenomicRanges-method}
\alias{ncol,GenomicRanges-method}
\alias{colnames,GenomicRanges-method}
\alias{colnames<-,GenomicRanges-method}

% Sequence methods:
\alias{[,GenomicRanges-method}
\alias{c,GenomicRanges-method}
\alias{length,GenomicRanges-method}
\alias{rev,GenomicRanges-method}
\alias{seqselect,GenomicRanges-method}
\alias{seqselect<-,GenomicRanges-method}
\alias{split,GenomicRanges-method}
\alias{window,GenomicRanges-method}

% show method:
\alias{show,GenomicRanges-method}


\title{GenomicRanges objects}

\description{
  The GenomicRanges class is a container for genomic locations and
  their associated annotations.
}

\details{
  The GenomicRanges class stores sequences of genomic locations and
  associated annotations. Each element in the sequence is comprised of a
  sequence name, an interval, a \link{strand}, and optional data values
  (e.g. score, GC content, etc.). This information is stored in four
  slots:
  \describe{
    \item{\code{seqnames}}{a 'character' \link[IRanges]{Rle} object
      containing the sequence names.}
    \item{\code{ranges}}{an \link[IRanges]{IRanges} object containing
      the ranges.}
    \item{\code{strand}}{a 'factor' \link[IRanges]{Rle} object containing
      the \link{strand} information.}
    \item{\code{values}}{a \link[IRanges]{DataFrame} object containing
      the annotation columns.}
  }
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicRanges(seqnames = Rle(), ranges = IRanges(),
                          strand = Rle(NA_character_, length(seqnames)), ...)}:
      Creates a GenomicRanges object.
      \describe{
        \item{\code{seqnames}}{Rle object or character vector containing the
          sequence names.}
        \item{\code{ranges}}{IRanges object containing the ranges.}
        \item{\code{strand}}{Rle object, character vector, or factor
          containing the strand information.}
        \item{\code{\ldots}}{Optional annotation columns for the \code{values}
          slot.}
      }
    }
  }
}

\section{Coercion}{
  In the code snippets below, \code{x} is a GenomicRanges object.

  \describe{
    \item{}{
      \code{as(from, "GenomicRanges")}:
      Creates a GenomicRanges object from a RangedData object.
    }
    \item{}{
      \code{as.data.frame(x, row.names = NULL, optional = FALSE)}:
      Creates a data.frame with columns \code{seqnames} (character), 
      \code{start} (integer), \code{end} (integer), \code{width} (integer),
      \code{strand} (factor), as well as the additional columns stored in
      \code{values(x)}.
    }
  }
}

\section{Accessors}{
  In the following code snippets, \code{x} is a GenomicRanges object.

  \describe{
    \item{}{
      \code{seqnames(x)}, \code{seqnames(x) <- value}:
      Get or set the sequence names. \code{value} can be an Rle object or
      character vector.
    }
    \item{}{
      \code{ranges(x)}, \code{ranges(x) <- value}:
      Get or set the ranges. \code{value} can be a Ranges object.
    }
    \item{}{
      \code{strand(x)}, \code{strand(x) <- value}:
      Get or set the strand. \code{value} can be an Rle object, character
      vector, or factor.
    }
    \item{}{
      \code{values(x)}, \code{values(x) <- value}:
      Get or set the optional data columns. \code{value} can be a DataFrame,
      data.frame object, or NULL.
    }
    \item{}{
      \code{names(x)}, \code{names(x) <- value}:
      Get or set the names of the elements.
    }
    \item{}{
      \code{length(x)}:
      Get the number of elements.
    }
  }
}

\section{Ranges methods}{
  In the following code snippets, \code{x} is a GenomicRanges object.

  \describe{
    \item{}{
      \code{start(x)}, \code{start(x) <- value}:
      Get or set \code{start(ranges(x))}.
    }
    \item{}{
      \code{end(x)}, \code{end(x) <- value}:
      Get or set \code{end(ranges(x))}.
    }
    \item{}{
      \code{width(x)}, \code{width(x) <- value}:
      Get or set \code{width(ranges(x))}.
    }
  }
}

\section{DataTable methods}{
  In the following code snippets, \code{x} is a GenomicRanges object.

  \describe{
    \item{}{
      \code{ncol(x)}:
      Get the number of column in \code{values(x)}.
    }
    \item{}{
      \code{colnames(x, do.NULL = TRUE, prefix = "col")},
      \code{colnames(x) <- value}:
      Get or set the column names of \code{values(x)}.
    }
  }
}

\section{Splitting and Combining}{
  In the code snippets below, \code{x} is a GenomicRanges object.

  \describe{
    \item{}{
      \code{append(x, values, after = length(x))}:
      Insert the \code{values} into \code{x} at the position given by
      \code{after}, where \code{x} and \code{values} are of the same
      class.
    }
    \item{}{
      \code{c(x, ...)}:
      Combine \code{x} and the GenomicRanges objects in \code{...} together.
      Any object in \code{...} must belong to the same class as \code{x},
      or to one of its subclasses, or must be \code{NULL}.
      The result is an object of the same class as \code{x}.
    }
    \item{}{
      \code{split(x, f, drop = FALSE)}:
      Splits \code{x} into a \link{GenomicRangesList}, according to \code{f},
      dropping elements corresponding to unrepresented levels if \code{drop}
      is \code{TRUE}.
    }
  }
}

\section{Subsetting}{
  In the code snippets below, \code{x} is a GenomicRanges object.

  \describe{
    \item{}{
      \code{x[i, j]}, \code{x[i, j] <- value}:
      Get or set elements \code{i} with optional values columns
      \code{values(x)[,j]}, where \code{i} can be missing; an NA-free
      logical, numeric, or character vector; or a 'logical' Rle object.
    }
    \item{}{
      \code{x[i,j] <- value}:
      Replace elements \code{i} and optional values columns \code{j} with
      \code{values}.
    }
    \item{}{
      \code{head(x, n = 6L)}:
      If \code{n} is non-negative, returns the first n elements of the
      GenomicRanges object.
      If \code{n} is negative, returns all but the last \code{abs(n)} elements
      of the GenomicRanges object.
    }
    \item{}{
      \code{rep(x, times, length.out, each)}:
      Repeats the values in \code{x} through one of the following conventions:
      \describe{
        \item{\code{times}}{Vector giving the number of times to repeat each
          element if of length \code{length(x)}, or to repeat the whole vector
          if of length 1.}
        \item{\code{length.out}}{Non-negative integer. The desired length of
          the output vector.}
        \item{\code{each}}{Non-negative integer.  Each element of \code{x} is
          repeated \code{each} times.}
      }
    }
    \item{}{
      \code{rev(x)}:
      Return a new object of the same class as \code{x} made of the original
      elements in the reverse order.
    }
    \item{}{
      \code{seqselect(x, start=NULL, end=NULL, width=NULL)}:
      Similar to \code{window}, except that multiple consecutive subsequences
      can be requested for concatenation. As such two of the three \code{start},
      \code{end}, and \code{width} arguments can be used to specify the
      consecutive subsequences. Alternatively, \code{start} can take a Ranges
      object or something that can be converted to a Ranges object like an
      integer vector, logical vector or logical Rle. If the concatenation of
      the consecutive subsequences is undesirable, consider using
      \code{\link{Views}}.
    }
    \item{}{
      \code{seqselect(x, start=NULL, end=NULL, width=NULL) <- value}:
      Similar to \code{window<-}, except that multiple consecutive subsequences
      can be replaced by a \code{value} whose length is a divisor of the number
      of elements it is replacing. As such two of the three \code{start},
      \code{end}, and \code{width} arguments can be used to specify the
      consecutive subsequences. Alternatively, \code{start} can take a Ranges
      object or something that can be converted to a Ranges object like an
      integer vector, logical vector or logical Rle.
    }
    \item{}{
      \code{subset(x, subset)}:
      Return a new object of the same class as \code{x} made of the subset
      using logical vector \code{subset}, where missing values are taken as
      \code{FALSE}.
    }
    \item{}{
      \code{tail(x, n = 6L)}:
      If \code{n} is non-negative, returns the last n elements of the
      GenomicRanges object.
      If \code{n} is negative, returns all but the first \code{abs(n)} elements
      of the GenomicRanges object.
    }
    \item{}{
      \code{window(x, start = NA, end = NA, width = NA, frequency = NULL, delta = NULL, ...)}:
      Extract the subsequence window from the GenomicRanges object using:
      \describe{
        \item{\code{start}, \code{end}, \code{width}}{The start, end, or width
          of the window. Two of the three are required.}
        \item{\code{frequency}, \code{delta}}{Optional arguments that specify
          the sampling frequency and increment within the window.}
      }
      In general, this is more efficient than using \code{"["} operator.
    }
    \item{}{
      \code{window(x, start = NA, end = NA, width = NA, keepLength = TRUE) <- value}:
      Replace the subsequence window specified on the left (i.e. the subsequence
      in \code{x} specified by \code{start}, \code{end} and \code{width})
      by \code{value}.
      \code{value} must either be of class \code{class(x)}, belong to a subclass
      of \code{class(x)}, be coercible to \code{class(x)}, or be \code{NULL}.
      If \code{keepLength} is \code{TRUE}, the elements of \code{value} are
      repeated to create a GenomicRanges object with the same number of elements
      as the width of the subsequence window it is replacing.
      If \code{keepLength} is \code{FALSE}, this replacement method can modify
      the length of \code{x}, depending on how the length of the left
      subsequence window compares to the length of \code{value}.
    }
  }
}

\author{P. Aboyoun}

\seealso{
  \link{GenomicRangesList},
  \link[IRanges]{Sequence},
  \link[IRanges]{Ranges},
  \link[IRanges]{Rle},
  \link[IRanges]{DataFrame}
}

\examples{
}
