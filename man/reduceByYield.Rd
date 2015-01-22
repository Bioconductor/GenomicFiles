\name{reduceByYield}
\alias{reduceByYield}

\title{
  Iterate through a BAM (or other) file, reducing output to a single result.
}
\description{

  Rsamtools files can be created with a \sQuote{yieldSize} argument that
  influences the number of records (chunk size) input at one time (see,
  e.g,. \code{\link[Rsamtools]{BamFile}}). \code{reduceByYield} iterates 
  through the file, processing each chunk and reducing it with previously 
  input chunks. This is a memory efficient way to process large data files,
  especially when the final result fits in memory.

}

\usage{
reduceByYield(X, YIELD,  MAP = identity, REDUCE = `+`, 
              DONE = function(x) is.null(x) || length(x) == 0L, 
              ..., parallel = FALSE, iterate = TRUE, init)
}

\arguments{

  \item{X}{A \code{\link[Rsamtools]{BamFile}} instance (or other class 
    for which \code{isOpen}, \code{open}, \code{close} methods are defined, 
    and which support extraction of sequential chunks).}

  \item{YIELD}{A function name or user-supplied function that operates 
    on \code{X} to produce a \code{VALUE} that is passed to \code{DONE} 
    and \code{MAP}. Generally \code{YIELD} will be a data extractor such as
    \code{readGAlignments}, \code{scanBam}, \code{yield}, etc. and
    \code{VALUE} is a chunk of data.
    \itemize{
      \item YIELD(X)
    }}

  \item{MAP}{A function of one or more arguments that operates on
    the chunk of data from \code{YIELD}.
    \itemize{
      \item MAP(VALUE, ...)
    }}

  \item{REDUCE}{A function of one (\code{iterate=FALSE} or two
    (\code{iterate=TRUE}) arguments, returning the reduction (e.g.,
    addition) of the argument(s). If missing, \code{REDUCE} is \code{c}
    (when \code{iterate=TRUE}) or \code{identity} when
    (when \code{iterate=FALSE}).
    \itemize{
      \item REDUCE(all.mapped, ...) ## iterate=FALSE
      \item REDUCE(x, y, ...)       ## iterate=TRUE
    }}

  \item{DONE}{A function of one argument, the \code{VALUE} output of 
    the most recent call to \code{YIELD(X, ...)}. If missing, \code{DONE} 
    is \code{function(VALUE) length(VALUE) == 0}.}

  \item{\dots}{Additional arguments, passed to \code{MAP}.}

  \item{iterate}{logical(1) determines whether the call to
    \code{REDUCE} is iterative (\code{iterate=TRUE}) or cumulative
    (\code{iterate=FALSE}).}

  \item{parallel}{logical(1) determines if the \code{MAP} step 
    is run in parallel. \code{bpiterate} is used under the hood
    and is currently supported for Unix/Mac only. For Windows machines, 
    \code{parallel} is ignored.}

  \item{init}{(Optional) Initial value used for \code{REDUCE} when
    \code{iterate=TRUE}.}

}

\details{
 
  When \code{iterate=TRUE}, \code{REDUCE} requires 2 arguments and is
  invoked with \code{init} and the output from the first call to
  \code{MAP}. If \code{init} is missing, it operates on the first two
  outputs from \code{MAP}.

  When \code{iterate=FALSE}, \code{REDUCE} requires 1 argument and is
  is invoked with a list containing a list containing all results from
  \code{MAP}.

}
\value{

  The value returned by the final invocation of \code{REDUCE}, or \code{init} 
  if provided and no data were yield'ed, or \code{list()} if \code{init} is 
  missing and no data were yield'ed.

}

\author{Martin Morgan \url{mtmorgan@fhcrc.org} and
        Valerie Obenchain \url{vobencha@fhcrc.org}}

\seealso{
  \itemize{
    \item \code{\link[Rsamtools]{BamFile}} and
          \code{\link[Rsamtools]{TabixFile}} for examples of `X`.
    \item \code{reduceByFile} and \code{reduceByRange}
  }
}
\examples{

if (all(require(RNAseqData.HNRNPC.bam.chr14) &&
        require(GenomicAlignments))) { 

  ## -----------------------------------------------------------------------
  ## Nucleotide frequency of mapped reads
  ## -----------------------------------------------------------------------
 
  ## In this example nucleotide frequency of mapped reads is computed
  ## for a single file. The MAP step is run in parallel and REDUCE 
  ## is iterative.

  ## Create a BamFile and set a 'yieldSize'.
  fl <- system.file(package="Rsamtools", "extdata", "ex1.bam")
  bf <- BamFile(fl, yieldSize=500)

  ## Define 'YIELD', 'MAP' and 'REDUCE' functions.
  YIELD <- function(X, ...) {
      flag = scanBamFlag(isUnmappedQuery=FALSE)
      param = ScanBamParam(flag=flag, what="seq")
      scanBam(X, param=param, ...)[[1]][['seq']]
  }
  MAP <- function(value, ...) {
      library(Biostrings) ## for alphabetFrequency()
      alphabetFrequency(value, collapse=TRUE)
  }
  REDUCE <- `+`        # add successive alphabetFrequency matrices

  ## 'parallel=TRUE' runs the MAP step in parallel and is currently
  ## implemented for Unix/Mac only.
  register(MulticoreParam(3))
  reduceByYield(bf, YIELD, MAP, REDUCE, parallel=TRUE)
 
  ## -----------------------------------------------------------------------
  ## Coverage
  ## -----------------------------------------------------------------------
 
  ## If sufficient resources are available coverage can be computed
  ## across several large BAM files by combining reduceByYield() with
  ## bplapply().

  ## Create a BamFileList with a few sample files and a Snow cluster
  ## with the same number of workers as files.
  bfl <- BamFileList(RNAseqData.HNRNPC.bam.chr14_BAMFILES[1:3])
  bpparam <- SnowParam(length(bfl))

  ## 'FUN' is run on each worker. Because these are Snow workers each
  ## variable used in 'FUN' must be explicitly passed. (This is not the case
  ## when using Multicore.)
  FUN <- function(bf, YIELD, MAP, REDUCE, parallel, ...) {
    library(GenomicAlignments) ## for readGAlignments()
    library(GenomicFiles)      ## for reduceByYield()
    reduceByYield(bf, YIELD, MAP, REDUCE, parallel=parallel)
  }
 
  ## Passing parallel=FALSE to reduceByYield() runs the MAP step in serial on 
  ## each worker. In this example, parallel dispatch is at the file-level 
  ## only (bplapply()).
  YIELD <- `readGAlignments`
  MAP <- function(value, ...) coverage(value)[["chr14"]] 
  REDUCE <- `+`
  bplapply(bfl, FUN, YIELD=YIELD, MAP=MAP, REDUCE=REDUCE,
           parallel=FALSE, BPPARAM = bpparam) 
}
}

\keyword{manip}
