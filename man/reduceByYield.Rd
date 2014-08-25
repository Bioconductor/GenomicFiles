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
reduceByYield(X, YIELD, MAP, REDUCE, 
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
    is run in parallel.}

  \item{init}{(Optional) Initial value used for \code{REDUCE} when
    \code{iterate=TRUE}.}

}

\details{
  
  When \code{iterate=TRUE}, \code{REDUCE} is initially invoked with
  either the \code{init} value and the value of the first call to
  \code{MAP} or, if \code{init} is missing, the values of the first two
  calls to \code{MAP}.

  When \code{iterate=FALSE}, \code{REDUCE} is invoked with a list
  containing a list with as many elements as there were calls to
  \code{MAP}.  Each element the result of an invocation of \code{MAP}.

}
\value{

  The return value is the value returned by the final invocation of
  \code{REDUCE}, or \code{init} if provided and no data were yield'ed,
  or \code{list()} if \code{init} is missing and no data were yield'ed.

}

\author{Martin Morgan \url{mtmorgan@fhcrc.org}}

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
  
  fl <- system.file(package="Rsamtools", "extdata", "ex1.bam")
  bf <- BamFile(fl, yieldSize=500) ## typically, yieldSize=1e6
  
  param <- ScanBamParam(
      flag=scanBamFlag(isUnmappedQuery=FALSE),
      what="seq")
  YIELD <- function(X, ...) scanBam(X, param, ...)[[1]][['seq']]
  MAP <- function(value, ...) 
      alphabetFrequency(value, collapse=TRUE)
  REDUCE <- `+`        # add successive alphabetFrequency matrices 
  reduceByYield(bf, YIELD, MAP, REDUCE, param=param, parallel=TRUE)
  
  ## -----------------------------------------------------------------------
  ## Coverage
  ## -----------------------------------------------------------------------
  
  ## reduceByYield() can be applied to multiple files by combining it
  ## with bplapply().
  
  ## FUN will be run on each worker; it contains the necessary arguments 
  ## to reduceByYield() as well as a call to the function itself.
  ## reduceByYield() could also be run in parallel (parallel=TRUE) 
  ## but in this example it is not.
  FUN <- function(bf) {
    library(GenomicAlignments)
    library(GenomicFiles)
    YIELD <- `readGAlignments`
    MAP <- function(value, ...) coverage(value)[["chr14"]] 
    REDUCE <- `+`
    reduceByYield(bf, YIELD, MAP, REDUCE)
  }
  
  ## BAM files are distributed across Snow workers and each worker applies
  ## reduceByYield().
  bfl <- BamFileList(RNAseqData.HNRNPC.bam.chr14_BAMFILES[1:3])
  bplapply(bfl, FUN, BPPARAM = SnowParam(3)) 
}
}

\keyword{manip}
