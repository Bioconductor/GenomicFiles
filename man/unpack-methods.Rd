\name{unpack}

\alias{unpack}
\alias{unpack,list,GRangesList-method}
\alias{unpack,List,GRangesList-method}

\title{Un-packing results obtained with a pack()ed group of ranges}

\description{
  \code{unpack} is a generic function with a few methods in order to allow
  easy inversion of results obtained with ranges created with \code{pack(x)}.
}

\usage{
\S4method{unpack}{list,GRangesList}(flesh, skeleton, ...)
\S4method{unpack}{List,GRangesList}(flesh, skeleton, ...)
}

\arguments{
  \item{flesh}{
    A \code{List} object obtained from querying a file with \code{skeleton}.
  }
  \item{skeleton}{
    The pack()ed \code{GRangesList} object used to query a file. 
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
    The \code{unpack} method converts results obtained with
    a group of packed ranges back to original form.
}

\value{
  A unpacked form of \code{flesh}.
}

\seealso{
  \itemize{
    \item \code{\link{pack}} for packing ranges.
  }
}

\examples{
  fl <- system.file("extdata", "ex1.bam", package = "Rsamtools")
  gr <- GRanges(c(rep("seq2", 3), "seq1"), 
                IRanges(c(75, 1, 100, 1), width = 2))

  ## Ranges are packed by order within chromosome.
  gr
  pack(gr)

  ## When 'inter_range_len' is specified, ranges are packed in
  ## groups around the inter-range gap. See ?pack for details. 
  pk <- pack(gr, inter_range_len = 25)
  pk

  ## We send the packed groups of ranges to workers to
  ## compute coverage. The result is a list of length 3,
  ## the same length as the packed ranges.
  FUN <- function(ranges, fl, param) {
      bamWhich(param) <- ranges
      coverage(fl, param=param)[ranges]
  } 
  dat <- bplapply(as.list(pk), FUN, fl = fl, param = ScanBamParam())
  length(dat)

  ## Each list element contains an RleList of the coverage for
  ## the corresponding packed ranges.
  dat[[1]]
  pk[[1]]

  ## unpack() returns a RleList of results in the same order as 
  ## the unpacked 'gr' object.
  unpack(dat, pk)
}

\keyword{methods}
