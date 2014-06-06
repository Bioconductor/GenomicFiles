\name{reduceBy}

\alias{reduceBy}
\alias{reduceByFile}
\alias{reduceByFile,GRanges-method}
\alias{reduceByFile,GRangesList-method}
\alias{reduceByRange}
\alias{reduceByRange,GRanges-method}
\alias{reduceByRange,GRangesList-method}

\title{Parallel computations across files or ranges}

\description{
  Computations are distributed across files or ranges 
  with the option to iteratively combine results.
}

\usage{
\S4method{reduceByFile}{GRanges}(ranges, files, MAPPER, 
    REDUCER, iterate=FALSE, ..., init)
\S4method{reduceByFile}{GRangesList}(ranges, files, MAPPER, 
    REDUCER, iterate=FALSE, ..., init)

\S4method{reduceByRange}{GRanges}(ranges, files, MAPPER, 
    REDUCER, iterate=FALSE, ..., init)
\S4method{reduceByRange}{GRangesList}(ranges, files, MAPPER, 
    REDUCER, iterate=FALSE, ..., init)
}

\arguments{
  \item{ranges}{
    A \code{GRanges} or \code{GrangesList} object.
  }
  \item{files}{
    A \code{character} vector or \code{List} of filenames.
  }
  \item{MAPPER}{
    A function executed on each worker. The signature must 
    contain two arguments; the first represents the range(s) 
    and the second the file(s). There is no restriction on 
    the argument names and additional arguments may be
    provided.

    \code{MAPPER = function(range, file, ...)}
  }
  \item{REDUCER}{
    An optional function that combines (reduces) output from the
    \code{MAPPER}. The reduction is performed in the distributed 
    step and combines data from a single worker. The value of 
    \code{iterate} affects when the \code{REDUCER} is applied and 
    the expected function signature. See details below.

    If no \code{REDUCER} is provided no reduction is performed.
    Results are the output directly from the \code{MAPPER}.

    \itemize{
      \item{reduce with iteration}{
        When \code{iterate=TRUE} the first argument to
        \code{REDUCER} is the accumulation of past results 
        and the second argument is the current result. There 
        is no restriction on argument names and additional 
        arguments may be supplied.

        \code{REDUCER = function(last, current, ...)}

        The \code{REDUCER} is applied after the \code{MAPPER}
        has processed each record or element.
      }
      \item{reduce without iteration}{
        When \code{iterate=FALSE} the first argument of
        \code{REDUCER} is the list of all output from
        the \code{MAPPER}. There is no restriction on the
        argument name and additional arguments may be supplied.

        \code{REDUCER = function(mapped, ...)}

        The \code{REDUCER} is applied after the \code{MAPPER}
        has processed all records or elements.
      }
    }
  }
  \item{iterate}{
    A logical indicating if the \code{REDUCER} function
    should be applied iteratively to the output of 
    \code{MAPPER}.

    Collapsing results iteratively is useful when the number of
    records to be procssed is large (maybe complete files) but
    the end result is a much reduced representation of all records.
    Iteratively applying \code{REDUCER} reduces the amount of
    data on each worker at any one time and can substantially
    reduce the memory footprint.
  }
  \item{init}{
    An (optional) initial value for \code{REDUCER}. Applicable when
    \code{iterate=TRUE}. \code{init} must be an object of the same 
    type as the elements returned from \code{MAPPER}. \code{REDUCER} 
    logically adds \code{init} to the start (when proceeding left 
    to right) or end of results obtained with \code{MAPPER}. 
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
  \code{reduceByFile} enables the extraction, manipulation and
  combination of ranges from a file. Data are combined within a
  file not across files. 

  \code{reduceByRanges} enables the extraction, manipulation and
  combination of a group of ranges across files. Data are combined
  across files not within a file.

  Both \code{MAPPER} and \code{REDUCER} are applied in a distributed 
  step. Currently there is no 'built-in' ability to combine results 
  from different workers.
}

\value{
  A \link{list} object.
}

\seealso{
  \itemize{
    \item \link{GenomicFileViews-class}
  }
}

\author{
  Martin Morgan <mtmorgan@fhcrc.org> and 
  Valerie Obenchain <vobencha@fhcrc.org>
}

\examples{

library(RNAseqData.HNRNPC.bam.chr14)
fls <- RNAseqData.HNRNPC.bam.chr14_BAMFILES  ## 8 bam files

## -----------------------------------------------------------------------
## Basics of reduceByFile() and reduceByRange():
## -----------------------------------------------------------------------

## This first example uses a MAPPER function only (no REDUCER). 

## Ranges of interest.
gr <- GRanges("chr14", IRanges(c(19100000, 106000000), width=1e7))

## The MAPPER function counts the number of junctions in each range
## (i.e., 'N' operations in the CIGAR).
MAPPER <- function(range, file, ...) {
    library(GenomicAlignments)
    param = ScanBamParam(which=range)
    gal = readGAlignments(file, param=param)
              table(njunc(gal))
} 

## reduceByFile() works within each file.
rbf <- reduceByFile(gr, fls, MAPPER)

## Output length matches the number of files and elementLengths 
## correspond to the number of ranges.
length(rbf)          ## 8 files
elementLengths(rbf)  ## 2 ranges

## Each list element contains a table of counts, one for each range.
rbf[[1]]

## In contrast, reduceByRange() extracts data across files.
rbr <- reduceByRange(gr, fls, MAPPER)

## Output length corresponds to the number of ranges.
length(rbr)          ## 2 ranges
elementLengths(rbr)  ## 8 files

## Each list element contains a table of counts, one for each file.
do.call(rbind, rbr[[1]])

## -----------------------------------------------------------------------
## Using a REDUCER:
## -----------------------------------------------------------------------

## The REDUCER function operates on the output of the MAPPER. It
## can be applied iteratively (after each step of the MAPPER) or 
## non-iteratively (after MAPPER is done). 

## In this example we compute coverage for a group of ranges across
## all files. The output of a call to coverage() is an RleList
## one is produced each time MAPPER is called and these RleLists
## accumulate on the worker. If the REDUCER is called 
## iteratively the 'current' result is collapsed with the 'last'
## resulting in a maximum of 2 RleLists on a worker at a time.

## Regions of interest.
gr <- GRanges("chr14", IRanges(c(62262735, 63121531, 63980327),
              width=214700))

## The MAPPER computes the pileups ...
MAPPER <- function(range, file, ...) {
    library(GenomicRanges)
    param = ScanBamParam(which=range)
    coverage(file, param=param)
} 

## and the REDUCER adds the results. When the REDUCER is applied
## iteratively the first 2 function arguments correspond to the
## 'last' and 'current' results from MAPPER.
REDUCER1 <- function(last, current, ...)
    c(last, current)

## The goal is to sum coverage across files we call reduceByRange().
cov1 <- reduceByFile(gr, fls, MAPPER, REDUCER1, iterate=TRUE)
cov1[[1]]

## If memory use is not a concern the REDUCER can be applied
## non-iteratively. A non-iterative REDUCER must take the
## MAPPED output as the first argument.
REDUCER2 <- function(mapped, ...)
    do.call(c, mapped)

## Call reduceByFile() with 'iterate=FALSE'.
cov2 <- reduceByFile(gr, fls, MAPPER, REDUCER2, iterate=FALSE)

## Results are the same as the with the iterative REDUCER.
cov2[[1]]

## -----------------------------------------------------------------------
## mean coverage per-file:
## -----------------------------------------------------------------------

## FIXME: better BigWig data
fl <- system.file("tests", "test.bw", package = "rtracklayer")
fls <- c(fl, fl, fl)
gr <- GRanges(c("chr2", "chr19", "chr19"),
              IRanges(c(1, 1400, 1700), width = 1000))

## MAPPER returns total coverage and counts (width) for each range.
MAPPER = function(range, file, ...) {
    rle <- import(file, selection=range, as="RleList")
    dat <- Filter(length, Views(rle, as(range, "RangesList")))

    c(n = as.numeric(sum(unlist(width(dat)))),
      sum = as.numeric(sum(sapply(dat, sum))))
}

## REDUCER sums coverage and counts for all ranges.
REDUCER = function(X, Y, ...)
    c(n = X[["n"]] + Y[["n"]], sum = X[["sum"]] + Y[["sum"]])

res <- reduceByFile(gr, fls, MAPPER, REDUCER, iterate=TRUE)
}

\keyword{methods}
