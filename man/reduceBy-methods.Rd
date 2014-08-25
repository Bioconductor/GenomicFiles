\name{reduceBy}

\alias{reduceBy}
\alias{reduceByFile}
\alias{reduceByFile,GRanges,ANY-method}
\alias{reduceByFile,GRangesList,ANY-method}
\alias{reduceByFile,GenomicFiles,missing-method}
\alias{reduceByRange}
\alias{reduceByRange,GRanges,ANY-method}
\alias{reduceByRange,GRangesList,ANY-method}
\alias{reduceByRange,GenomicFiles,missing-method}

\title{Parallel computations across files or ranges}

\description{
  Computations are distributed across files or ranges 
  with the option to iteratively combine results.
}

\usage{
\S4method{reduceByFile}{GRanges,ANY}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)
\S4method{reduceByFile}{GRangesList,ANY}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)
\S4method{reduceByFile}{GenomicFiles,missing}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)

\S4method{reduceByRange}{GRanges,ANY}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)
\S4method{reduceByRange}{GRangesList,ANY}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)
\S4method{reduceByRange}{GenomicFiles,missing}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)
}

\arguments{
  \item{ranges}{
    A \code{GRanges}, \code{GrangesList} or \code{GenomicFiles}
    object. When \code{ranges} is a \code{GenomicFiles}, the
    \code{files} argument is missing.
  }
  \item{files}{
    A \code{character} vector or \code{List} of filenames.
    When \code{ranges} is a \code{GenomicFiles} the \code{files}
    argument is missing; both ranges and files are extracted 
    from the object.
  }
  \item{MAP}{
    A function executed on each worker. The signature must 
    contain two arguments; the first represents the range(s) 
    and the second the file(s). There is no restriction on 
    the argument names and additional arguments may be
    provided.

    \itemize{
      \item \code{MAP = function(range, file, ...)}
    }
  }
  \item{REDUCE}{
    An optional function that combines (reduces) output from the
    \code{MAP}. The first argument is the list output from MAP,
    additional arguments may be supplied. There are no restrictions
    on argument names.

    \itemize{
      \item \code{REDUCE = function(mapped, ...)}
    }

    Reduction combines data from a single worker and is always
    performed as part of the distributed step. When \code{iterate=TRUE} 
    \code{REDUCE} is applied after each \code{MAP} step;
    depending on the nature of \code{REDUCE}, iterative reduction
    can substantially decrease the data stored in memory. When 
    \code{iterate=FALSE} reduction is applied to the list of \code{MAP} 
    output applied to all files / ranges.

    When \code{REDUCE} is missing, output is a list from \code{MAP}.
  }
  \item{iterate}{
    A logical indicating if the \code{REDUCE} function
    should be applied iteratively to the output of 
    \code{MAP}. When \code{REDUCE} is missing \code{iterate}
    is set to FALSE. 

    Collapsing results iteratively is useful when the number of
    records to be procssed is large (maybe complete files) but
    the end result is a much reduced representation of all records.
    Iteratively applying \code{REDUCE} reduces the amount of
    data on each worker at any one time and can substantially
    reduce the memory footprint.
  }
  \item{summarize}{
    A logical indicating if results should be returned as a
    \code{SummarizedExperiment} object instead of a list.
    \code{SummarizedExperiment} requires matching dimensions
    across rows and columns of the slots. Becasue a \code{REDUCE} 
    collapses one dimension (either ranges or files) the result
    cannnot be put in a \code{SummarizedExperiment}. When a 
    \code{REDUCE} is provided \code{summarize} is ignored
    (i.e., set to FALSE).
  }
  \item{init}{
    An (optional) initial value for \code{REDUCE} when
    \code{iterate=TRUE}. \code{init} must be an object of the same 
    type as the elements returned from \code{MAP}. \code{REDUCE} 
    logically adds \code{init} to the start (when proceeding left 
    to right) or end of results obtained with \code{MAP}. 
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\details{
  The reduceBy* functions offer two approaches to working with
  data subsets from multiple files. \code{reduceByFile} enables the 
  extraction, manipulation and combination of data within files
  while \code{reduceByRanges} works across files.

  Both \code{MAP} and \code{REDUCE} functions can be provided
  but only the \code{MAP} is required. The first two arguments to
  \code{MAP} are `range' and `file'; the first argument to 
  \code{REDUCE} is the list of output from the \code{MAP}.

  Both \code{MAP} and \code{REDUCE} are applied in a distributed 
  step. Currently there is no 'built-in' ability to combine results 
  across workers in the distributed step.
}

\value{
  Output is a \code{list} when \code{summarize=FALSE} (default)
  and a \code{SummarizedExperiment} when \code{summarize=TRUE}.
  Note that if \code{REDUCE} is provided \code{summarize} is
  ignored (i.e., set to FALSE).

  When \code{ranges} is a \code{GenomicFiles} object and 
  \code{summarize=TRUE}, data from \code{rowData}, \code{colData} and 
  \code{exptData} are transferred to the \code{SummarizedExperiment}.
}

\seealso{
  \itemize{
    \item \link{GenomicFiles-class}
  }
}

\author{
  Martin Morgan <mtmorgan@fhcrc.org> and 
  Valerie Obenchain <vobencha@fhcrc.org>
}

\examples{

if (all(require(RNAseqData.HNRNPC.bam.chr14) &&
        require(GenomicAlignments))) {
  fls <- RNAseqData.HNRNPC.bam.chr14_BAMFILES  ## 8 bam files
  
  ## -----------------------------------------------------------------------
  ## Basics of reduceByFile() and reduceByRange():
  ## -----------------------------------------------------------------------
  
  ## In this first example we provide a MAP only (no REDUCE).
  
  ## Ranges of interest.
  gr <- GRanges("chr14", IRanges(c(19100000, 106000000), width=1e7))
  
  ## The MAP counts the number of junctions in each range
  ## (i.e., 'N' operations in the CIGAR).
  MAP <- function(range, file, ...) {
      library(GenomicAlignments)
      param = ScanBamParam(which=range)
      gal = readGAlignments(file, param=param)
                table(njunc(gal))
  } 
  
  ## Length of the output corresponds to the number of files and 
  ## the elementLengths to the number of ranges.
  rbf <- reduceByFile(gr, fls, MAP)
  length(rbf)          ## 8 files
  elementLengths(rbf)  ## 2 ranges
  
  ## Each list element contains a table of counts, one for each range.
  rbf[[1]]
  
  ## In contrast, reduceByRange() extracts data across files.
  rbr <- reduceByRange(gr, fls, MAP)
  
  ## Output length corresponds to the number of ranges.
  length(rbr)          ## 2 ranges
  elementLengths(rbr)  ## 8 files
  
  ## Each list element contains a table of counts, one for each file.
  do.call(rbind, rbr[[1]])
  
  ## Output a SummarizedExperiment instead of list:
  se <- reduceByRange(gr, fls, MAP, summarize=TRUE)
  assays(se)
  
  ## -----------------------------------------------------------------------
  ## Computing coverage across files:
  ## -----------------------------------------------------------------------
  
  ## Use reduceByRange() to compute coverage for a group of ranges
  ## across files.
  
  ## Regions of interest.
  gr <- GRanges("chr14", IRanges(c(62262735, 63121531, 63980327),
                width=214700))
  
  ## The MAP computes the pileups ...
  MAP <- function(range, file, ...) {
      library(GenomicRanges)
      param = ScanBamParam(which=range)
      coverage(file, param=param)[range]
  } 
  
  ## and the REDUCE adds the last and current results. 
  REDUCE <- function(mapped, ...)
      Reduce("+", mapped)
  
  ## Each call to coverage() produces an RleList which accumulate 
  ## on the workers. When the REDUCE is applied iteratively the
  ## 'current' result is collapsed with the 'last' resulting in a 
  ## maximum of 2 RleLists on a worker at a time.
  cov1 <- reduceByRange(gr, fls, MAP, REDUCE, iterate=TRUE)
  cov1[[1]]
  
  ## If memory use is not a concern (or if MAP output
  ## is not large) the REDUCE can be applied non-iteratively. 
  cov2 <- reduceByRange(gr, fls, MAP, REDUCE, iterate=FALSE)
  
  ## Results match those obtained with the iterative REDUCE.
  cov2[[1]]
  
  ## -----------------------------------------------------------------------
  ## Organizing runs with the GenomicFiles class:
  ## -----------------------------------------------------------------------
  
  ## The GenomicFiles class is a light-weight form of SummarizedExperiment
  ## that does not have an 'assays' slot. 
  colData <- DataFrame(method=rep("RNASeq", length(fls)),
                       format=rep("bam", length(fls)))
  gf <- GenomicFiles(files=fls, rowData=gr, colData=colData)
  gf
  
  ## The object can be subset on ranges or files for different
  ## experimental runs.
  dim(gf)
  gf_sub <- gf[2, 3:4]
  dim(gf_sub)
  
  ## When summarize = TRUE and no REDUCE is provided the reduceBy* 
  ## functions output a SummarizedExperiment object.
  se <- reduceByFile(gf, MAP=MAP, summarize=TRUE)
  se
  
  ## Data from the rowData, colData and exptData slots in the
  ## GenomicFiles are transferred to the SummarizedExperiment.
  colData(se)
  
  ## Results are in the assays slot.
  assays(se) 
}
}

\keyword{methods}
