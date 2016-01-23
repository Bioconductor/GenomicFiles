\name{unpack}

\alias{unpack}
\alias{unpack,list,GRangesList-method}
\alias{unpack,List,GRangesList-method}

\title{Un-pack results obtained with a pack()ed group of ranges}

\description{
  \code{unpack} returns results obtained with pack()ed ranges
  to the geometry of the original, unpacked ranges.
}

\usage{
\S4method{unpack}{list,GRangesList}(flesh, skeleton, ...)
\S4method{unpack}{List,GRangesList}(flesh, skeleton, ...)
}

\arguments{
  \item{flesh}{
    A \code{List} object to be unpacked; the result from
    querying a file with \code{skeleton}.
  }
  \item{skeleton}{
    The \code{GRangesList} created with `pack(x)`. 
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
    \code{unpack} returns a \code{List} obtained with
    packed ranges to the geometry and order of the original,
    unpacked ranges. 
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

  ## Ranges are packed by order within chromosome and grouped
  ## around gaps greater than 'inter_range_len'. See ?pack for details.
  pk <- pack(gr, inter_range_len = 25)

  ## FUN computes coverage for the range passed as 'rng'.
  FUN <- function(rng, fl, param) {
      requireNamespace("GenomicAlignments") ## for bamWhich() and coverage()
      Rsamtools::bamWhich(param) <- rng
      GenomicAlignments::coverage(Rsamtools::BamFile(fl), param=param)[rng]
  } 

  ## Compute coverage on the packed ranges.
  dat <- bplapply(as.list(pk), FUN, fl = fl, param = ScanBamParam())

  ## The result list contains RleLists of coverage.
  lapply(dat, class)

  ## unpack() transforms the results back to the order of 
  ## the original ranges (i.e., unpacked 'gr').
  unpack(dat, pk)
}

\keyword{methods}
