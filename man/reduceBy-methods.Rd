\name{reduceBy}

\alias{reduceBy}
\alias{reduceByFile}
\alias{reduceByFile,GenomicFileViews-method}
\alias{reduceByRange}
\alias{reduceByRange,GenomicFileViews-method}

\title{Distribute parallel computations across files}

\description{
  Distribute parallel computations across files with
  the option to provide MAP and REDUCE functions
  across ranges and / or files. 
}

\usage{
\S4method{reduceByFile}{GenomicFileViews}(X, MAP, REDUCE, ..., 
    init, ITERATE=FALSE)

\S4method{reduceByRange}{GenomicFileViews}(X, MAP, REDUCE, ..., 
    init, ITERATE=FALSE)
}

\arguments{
  \item{X}{
    A \code{GenomicFileViews} object.
  }
  \item{MAP}{
    A function applied to the files in \code{fileList} or ranges in 
    \code{fileRange}. \code{MAP} should have (at minimum) the two 
    arguments \code{RANGE} and \code{FILE}. Signature is as follows:

    \code{MAP = function(FILE, RANGE, ...)}

    \code{MAP} can be used without a \code{REDUCE} function.
    In this case results from \code{reduceByRange} and 
    \code{reduceByFile} are the same but returned in different order 
    with different geometry (lists and list elements are different 
    lengths).
  }
  \item{REDUCE}{
    A function applied to the output of \code{MAP} which (usually)
    performs an aggregation. The function signature of 
    \code{REDUCE} depends on when aggregation occurs. When
    \code{ITERATE} is FALSE \code{REDUCE} is applied after the mapping 
    step is complete; the signature requires a single argument
    (\code{MAPPED} below). When \code{ITERATE} is TRUE data are
    aggregated while \code{MAP} iterates through the files or ranges;
    the signature requires two arguments (\code{X} and \code{Y} below). 

    \itemize{
      \item{ITERATE = FALSE :}{
        MAPPED is the complete output from \code{MAP}.

        \code{REDUCE = function(MAPPED, ...)}
      }
      \item{ITERATE = TRUE :}{
        X is the result from the last reduce step and Y is the 
        most recent yield from \code{MAP}.

        \code{REDUCE = function(X, Y, ...)}
      }
    }
  }
  \item{init}{
    An (optional) initial value for \code{REDUCE}. Applicable when
    \code{ITERATE = TRUE}. \code{init} must be an object of the same 
    type as the elements returned from \code{MAP}. \code{REDUCE} 
    logically adds \code{init} to the start (when proceeding left 
    to right) or end of results obtained with \code{MAP}. 
    See ?\code{Reduce} for details. 
  }
  \item{ITERATE}{
    A \code{logical} specifying if reduction should be
    performed during the mapping step (TRUE) or after (FALSE).
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
  Computations are distributed in parallel by file or by range by
  specifying \code{reduceByFile} or \code{reduceByRange}. In the
  distributed step, both \code{reduceByFile} and \code{reduceByRange} 
  use MAP and REDUCE functions to further process and combine the data. 
  \code{MAP} and \code{REDUCE} are based on the list-processing combinators 
  from functional programming in the same spirit as \code{Map} and 
  \code{Reduce} in \code{base} R. The all-caps designation of \code{MAP} 
  and \code{REDUCE} is used to distinguish the functions in 
  \code{GenomicFileViews}.

  In the case of \code{reduceByRange} the mapping occurs between a 
  single range and all files. The reduce step combines the results
  across files within the single range. With \code{reduceByFile} the 
  mapping occurs between a single file and all ranges. The reduce step 
  combines results across ranges within a single file.

  The \code{ITERATE} argument controls when the \code{REDUCE}
  function is applied. When \code{TRUE}, \code{REDUCE} is applied 
  while \code{MAP} iterates through the files or ranges; when 
  \code{FALSE}, \code{REDUCE} is applied after \code{MAP} has 
  processed all files or ranges. 

  \itemize{
    \item{reduceByRange}{
      Computations are distributed in parallel over the ranges in
      \code{fileRange}. Each worker lapply()s the \code{MAP} function 
      to a single range over all files. The output is a list the 
      length of \code{fileRange} with each list element the length of 
      \code{fileList}.

      \code{ITERATE = FALSE}:
      Each worker applies \code{REDUCE} to \code{MAP} output from 
      a single range, across all files; independent combination 
      across files.

      \code{ITERATE = TRUE}:
      Each worker applies \code{REDUCE} as \code{MAP} iterates through 
      each file, single range; dependent combination across files.
    }
    \item{reduceByFile}{
      Computations are distributed in parallel by file. Each worker 
      lapply()s the \code{MAP} function over all ranges in a single 
      file. The output is a list the length of \code{fileList} with 
      each list element the length of \code{fileRange}.

      \code{ITERATE = FALSE}:
      Each worker applies \code{REDUCE} to \code{MAP} output from 
      a single file, across all ranges; independent combination 
      across ranges.

      \code{ITERATE = TRUE}:
      Each worker applies \code{REDUCE} as \code{MAP} iterates through 
      each range, a single file; dependent combination across ranges.
    }
  }
}

\value{
  A \link{list} object.
}

\seealso{
  \itemize{
    \item \link{GenomicFileViews-class}
    \item \link{Map}
    \item \link{Reduce}
  }
}

\author{
  Martin Morgan <mtmorgan@fhcrc.org> and 
  Valerie Obenchain <vobencha@fhcrc.org>
}

\examples{
if (.Platform$OS.type != "windows") {

  bw <- system.file("tests", "test.bw", package = "rtracklayer")
 
  ## -----------------------------------------------------------------------
  ## MAP only:
  ## -----------------------------------------------------------------------
  gr <- GRanges(c("chr2", "chr19", "chr19"),
                IRanges(c(1, 1400, 1700), width=10))
  bwfv <- BigWigFileViews(c(bw, bw), fileRange=gr)
  MAP = function(FILE, RANGE, ...) {
      import(FILE, selection=RANGE, as="NumericList")
  }

  ## reduceByRange(): 
  ## List length is the number of ranges, list element lengths are 
  ## the number of files.
  reduceByRange(bwfv, MAP)
 
  ## reduceByFile(): 
  ## List length is the number of files, list element lengths are 
  ## the number of ranges.
  reduceByFile(bwfv, MAP)

  ## -----------------------------------------------------------------------
  ## MAP, REDUCE and ITERATE:
  ## -----------------------------------------------------------------------
  ## This simple example demonstrates the use of ITERATE to
  ## control when output from MAP is aggregated.

  gr <- tileGenome(seqlengths(BigWigFile(bw)), 
                   tilewidth = 1e7,
                   cut.last.tile.in.chrom = TRUE)
  bwv <- BigWigFileViews(c(bw, bw), fileRange = gr)
  MAP = function(FILE, RANGE, ...) TRUE
 
  ## REDUCE is applied after MAP is complete.
  REDUCE = function(MAPPED, ...) sum(unlist(MAPPED))
  result1 <- reduceByRange(bwv, MAP, REDUCE)
 
  ## REDUCE is applied while MAP iterates though each file.
  REDUCE = function(X, Y, ...) X + Y
  result2 <- reduceByRange(bwv, MAP, REDUCE, ITERATE = TRUE)
 
  stopifnot(identical(result1, result2))
 
  ## -----------------------------------------------------------------------
  ## t-test at point locations across files:
  ## -----------------------------------------------------------------------
 
  library(genefilter)
  gr <- tileGenome(seqlengths(BigWigFile(bw)), 
                   tilewidth = 1e7,
                   cut.last.tile.in.chrom=TRUE)
  bwv <- BigWigFileViews(c(bw, bw, bw, bw), fileRange = gr)
  grp <- factor(c(1, 2, 2, 1))
 
  MAP = function(FILE, RANGE, ...) {
      stopifnot(length(RANGE) == 1)
      v <- Views(import(FILE, selection = RANGE, asRle = TRUE),
                 as(RANGE, "RangesList"))
      as.vector(Filter(length, v)[[1]][[1]])
  }
  REDUCE = function(MAPPED, ..., grp) {
      m <- simplify2array(MAPPED)
      idx <- which(rowSums(m) != 0)
      df <- rowttests(m[idx,], grp)
      cbind(offset=idx - 1, df)
  }
  result3 <- reduceByRange(bwv, MAP, REDUCE, grp = grp)
 
  ## Each list element contains the data.frame output
  ## from rowttests():
  result3[[1]]
 
  ## -----------------------------------------------------------------------
  ## mean per-file:
  ## -----------------------------------------------------------------------
 
  gr <- GRanges(c("chr2", "chr19", "chr19"),
                IRanges(c(1, 1400, 1700), width = 1000))
  bwv <- BigWigFileViews(c(bw, bw, bw), fileRange = gr)
 
  MAP = function(FILE, RANGE, ...) {
      v <- Filter(length, {
          Views(import(FILE, selection=RANGE, asRle = TRUE),
                as(RANGE, "RangesList"))
      })
      c(n = as.numeric(sum(unlist(width(v)))),
        sum = as.numeric(sum(sapply(v, sum))))
  }
  REDUCE = function(X, Y, ..., FILE)
      c(n = X[["n"]] + Y[["n"]], sum = X[["sum"]] + Y[["sum"]])
 
  result4 <- reduceByFile(bwv, MAP, REDUCE, ITERATE = TRUE)
}
}

\keyword{methods}
