\name{pack}

\alias{pack}
\alias{pack,GRanges-method}

\title{Range transformations of a GenomicRanges object for 
       optimal file queries
}

\description{
  Given a \code{GRanges} object, \code{pack} produces
  a like object obtained by grouping and/or concatenating ranges that
  meet specific criteria such as density, distance between ranges
  and max length. 
}

\usage{
\S4method{pack}{GRanges}(x, ..., 
    pack_type = c("all", "density", "distance", "length"),
    extraction_type = c("NumericList", "RleList"),
    max_len = 1e9, max_inter_range_len = 1e7, width = 1e7,
    density_ratio = 0.25) 
}

\arguments{
  \item{x}{
    A \code{GRanges} object.
  }
  \item{pack_type}{
    The character string specifying the pack method.
    Options are \code{density}, \code{distance}, \code{length}
    and \code{all} (default). See the `details' section for 
    descriptions of the packing methods.
  }
  \item{extraction_type}{
    The data format resulting from a query with the packed ranges.
    When \code{as=NumericList} the result is a NumericList the same 
    length as the packed ranges. For RleList, the result is an 
    RleList with one list element for each chromosome which is not
    the same length as the packed ranges. To extract the result
    in the same form as the original \code{x} use \code{\link{unpack}}.
  }
  \item{max_len}{
    A numeric specifying the max length allowed for ranges in \code{x}.
    Applicable when \code{pack_type=all} or \code{pack_type=length}.
    See Details section for packing description.
  }
  \item{max_inter_range_len}{
    A numeric specifying the max length allowed between ranges in \code{x}.
    Applicable when \code{pack_type=all} or \code{pack_type=distance}.
    See Details section for packing description.
  }
  \item{width}{
    A numeric specifying the width used in the \code{tile} function.
    Applicable when \code{pack_type=all} or \code{pack_type=length}.
  }
  \item{density_ratio}{
    A numeric between 0 and 1, specifying the ratio of coverage
    of the total range, for a set of ranges to be considered 'dense'.
    Applicable when \code{pack_type=all} or \code{pack_type=density}.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
    The \code{pack} method attempts to re-package ranges in optimal
    form for making file queries (i.e., data extraction). 

    \code{pack} checks one or all of the following qualities:
    \itemize{
      \item density: Ranges are dense if they cover more than
            \code{density_ratio} of the total range of \code{x}. Dense ranges
            are packed as a single max range.
 
      \item distance: Ranges are distant if the distance between any
            of the individual ranges exceeds \code{max_inter_range_len}.
            Distant ranges are packed in chunks (bins) around the gap
            separating the distant ranges.

      \item length: A range is long if it exceeds \code{max_len}.
            Long ranges are packed as sub-ranges using the
            \code{tile} function and \code{width} argument.
    } 
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
if (.Platform$OS.type != "windows") {
  gr1 <- GRanges("chr1", IRanges(1:5*5, width=3)) ## dense
  gr2 <- GRanges("chr2", IRanges(c(1:3, 30000:30003), width=1000)) ## distant
  gr3 <- GRanges("chr3", IRanges(c(1:2), width=c(45, 1e8))) ## long
 
  ## Dense ranges are converted to a single max range.
  p1 <- pack(gr1, pack_type="density")
  p1
 
  ## Ranges separated by > max_inter_range_len are grouped 
  ## into chunks defined by endpoints of the gap.
  p2 <- pack(gr2, pack_type="distance")
  p2
 
  ## Ranges exceeding 'max_len' are broken into chunks with
  ## the tile() function. 'width' defines the width of the
  ## sub-ranges created with tile().
  width(gr3)
  p3 <- pack(gr3, pack_type="length", max_len=1e7, width=1e6)
  p3
 
  ## When 'pack_type' is "all" (or not specified), ranges are checked 
  ## for all 3 characteristics, density, distance and length. 
  gr <- suppressWarnings(c(gr1, gr2, gr3))
  p_all <- pack(gr, pack_type="all", max_len=1e7, width=1e6)
  stopifnot(length(p_all) == sum(length(p1), length(p2), length(p3)))
}
}

\keyword{methods}
