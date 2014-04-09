\name{unpack}

\alias{unpack}
\alias{unpack,GRanges,RleList-method}
\alias{unpack,GRanges,NumericList-method}

\title{Un-packing results obtained with a pack()ed group of ranges}

\description{
  \code{unpack} is a generic function with a few methods in order to allow
  easy inversion of results obtained with ranges created with \code{pack(x)}.
}

\usage{
\S4method{unpack}{GRanges,NumericList}(x, y, ...)
\S4method{unpack}{GRanges,RleList}(x, y, ...)
}

\arguments{
  \item{x}{
    A \code{GRanges} object that was pack()ed and subsequently 
    used to query a file.
  }
  \item{y}{
    A \code{NumericList} or \code{RleList} obtained from querying
    a file with the pack()ed form of \code{x}.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
    The \code{unpack} method converts results obtained with
    a group of packed ranges back to original form. In
    other words, the data in \code{y} will be reassembled to
    match the geometry (ranges) in \code{x}.

    This function is under active development.
}

\value{
  A \code{NumericList} or \code{RleList} object.
}

\seealso{
  \itemize{
    \item \code{\link{pack}} for packing ranges.
  }
}

\examples{
if (.Platform$OS.type != "windows") {
  ## ---------------------------------------------------------------------
  ## Currently there are a limited number of objects (RleList and
  ## NumericList) that can be unpack()ed. Development this far has 
  ## been with BamFileView objects and is a work in progress.
  ## ---------------------------------------------------------------------

  fl <- system.file("tests", "test.bw", package = "rtracklayer")
  gr <- GRanges("chr2", IRanges(1:100, width=2))

  ## Dense ranges are packed as a single max range:
  pk <- pack(gr)

  ## Data comes back as a NumericList with the
  ## same geometry as the packed ranges: 
  dat <- import(fl, which=pk, as="NumericList")
  dat

  ## unlist() returns a NumericList in the same form as 'x': 
  unlist(gr, dat)
}
}

\keyword{methods}
