\name{pack}

\alias{pack}
\alias{isPacked}
\alias{pack,GRanges-method}

\title{Range transformations of a \code{GenomicRanges} object for 
       optimal file queries.
}

\description{
  Given a \code{GRanges} object, \code{pack} produces a \code{GRangesList}
  of the same ranges grouped and re-ordered.
}

\usage{
\S4method{pack}{GRanges}(x, ..., range_len = 1e9, inter_range_len = 1e7)
}

\arguments{
  \item{x}{
    A \code{GRanges} object.
  }
  \item{range_len}{
    A numeric specifying the max length allowed for ranges in \code{x}.
  }
  \item{inter_range_len}{
    A numeric specifying the max length allowed between ranges in \code{x}.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
  \subsection{Packing ranges}{
    The \code{pack} method attempts to re-package ranges in optimal
    form for extracting data from files. Ranges are not modified (made 
    shorter or longer) but re-ordered and / or re-grouped according 
    to the following criteria. 

    \itemize{
      \item order: Ranges are ordered by genomic position within
            chromosomes.
 
      \item distance: Ranges separted by a distance greater than
            the \code{inter_range_len} are packed in groups around 
            the gap separating the distant ranges.

      \item length: Ranges longer than \code{range_len} are packed 
            `individually' (i.e., retrived from the file as a single
             range vs grouped with other ranges).
    }
  }
  \subsection{Utilities}{
  \describe{
    \item{}{
      \code{isPacked(x, ...)}:
      Returns a logical indicating if the ranges in \code{x} are packed.
      \code{x} must be a \code{GRangesList} object.
    }
  }}
}

\value{
  A \code{GRanges} object.
}

\seealso{
  \itemize{
    \item \code{\link{unpack}} for unpacking the result obtained
    with `packed' ranges. 
  }
}

\examples{
## Ranges are ordered by position within chromosome.
gr1 <- GRanges("chr1", IRanges(5:1*5, width = 3)) 
pack(gr1)

## Ranges separated by > inter_range_len are partitioned
## into groups defined by the endpoints of the gap.
gr2 <- GRanges("chr2", IRanges(c(1:3, 30000:30003), width = 1000))
pack(gr2, inter_range_len = 20000)

## Ranges exceeding 'range_len' are isolated in a single element
## of the GRangesList.
gr3 <- GRanges("chr3", IRanges(c(1:4), width=c(45, 1e8, 45, 45)))
width(gr3)
pack(gr3, range_len = 1e7)
}

\keyword{methods}
