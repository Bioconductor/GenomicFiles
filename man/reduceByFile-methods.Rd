\name{reduceByFile}

\alias{reduceByFile}
\alias{reduceByFile,GRanges,ANY-method}
\alias{reduceByFile,GRangesList,ANY-method}
\alias{reduceByFile,GenomicFiles,missing-method}

\alias{reduceFiles}

\title{Parallel computations by files}

\description{
  Computations are distributed across files the with option to
  iteratively combine results.
}

\usage{
\S4method{reduceByFile}{GRanges,ANY}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)
\S4method{reduceByFile}{GRangesList,ANY}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)
\S4method{reduceByFile}{GenomicFiles,missing}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)

reduceFiles(ranges, files, MAP, REDUCE, ..., iterate=TRUE, init)
}

\arguments{
  \item{ranges}{
    A \code{GRanges}, \code{GrangesList} or \code{GenomicFiles}
    object.

    A \code{GRangesList} implies a grouping of the ranges; \code{MAP}
    is applied to each element of the \code{GRangesList} vs each range
    when \code{ranges} is a \code{GRanges}.

    When \code{ranges} is a \code{GenomicFiles} the \code{files}
    argument is missing; both ranges and files are extracted 
    from the object.
  }
  \item{files}{
    A \code{character} vector or \code{List} of filenames. A \code{List}
    implies a grouping of the files; \code{MAP} is applied to each
    element of the \code{List} vs each file individually.
  }
  \item{MAP}{
    A function executed on each worker. The signature must contain a
    minimum of two arguments representing the ranges and files. There is no 
    restriction on argument names and additional arguments can be provided.

    \itemize{
      \item \code{MAP = function(range, file, ...)}
    }
  }
  \item{REDUCE}{
    An optional function that combines output from the \code{MAP} step. The
    signature must contain at least one argument representing the list
    output from \code{MAP}. There is no restriction on argument names and 
    additional arguments can be provided.

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
    records to be processed is large (maybe complete files) but
    the end result is a much reduced representation of all records.
    Iteratively applying \code{REDUCE} reduces the amount of
    data on each worker at any one time and can substantially
    reduce the memory footprint.
  }
  \item{summarize}{
    A logical indicating if results should be returned as a
    \code{SummarizedExperiment} object instead of a list;
    data are returned in the \code{assays} slot named `data`. 
    This argument applies to \code{reduceByFile} only.

    When \code{REDUCE} is provided \code{summarize} is ignored
    (i.e., set to FALSE). A \code{SummarizedExperiment} requires the number 
    of rows in \code{rowData} and \code{assays} to match. Because \code{REDUCE} 
    collapses the data across ranges, the dimension of the result no longer
    matches that of the original ranges.
  }
  \item{init}{
    An optional initial value for \code{REDUCE} when
    \code{iterate=TRUE}. \code{init} must be an object of the same 
    type as the elements returned from \code{MAP}. \code{REDUCE} 
    logically adds \code{init} to the start (when proceeding left 
    to right) or end of results obtained with \code{MAP}. 
  }
  \item{\dots}{
    Arguments passed to other methods.
    \itemize{
      \item{nchunk}{
        Integer or NA used to control chunking of ranges in \code{reduceFiles}. 
        When length(ranges) is greater than 1e5, ranges are broken into chunks 
        with GenomicFiles:::.chunkIndex(). To override the internal algorithm
        set \code{nchunk} to desired integer length;to disable chunking 
        set \code{nchunk} to NA. For example, if length(ranges) is 1e6
        and \code{nchunk=3}, ranges are broken into 3 groups of approximately 
        equal length 333334, 333334, and 333332.
      }
    }
  }
}

\details{
  \code{reduceByFile} extracts, manipulates and combines multiple ranges 
  within a single file. Each file is sent to a worker where \code{MAP} is 
  invoked on each file / range combination. This approach allows multiple 
  ranges extracted from a single file to be kept separate or combined with 
  \code{REDUCE}.
  
  In contrast, \code{reduceFiles} does not iterate through the individual
  ranges but instead treats them as a group. \code{MAP} is invoked
  once for each file using all ranges as the \code{range} argument.
  In general, \code{REDUCE} does not play a significant role in 
  \code{reduceFiles} because \code{MAP} is only called once on each worker. 
  The exception is when the number of ranges exceeds 1e5 
  (i.e., length(ranges) > 1e5). In this case ranges are divided into chunks 
  and iterated through on each worker. To override the internal chunk 
  algorithm set \code{nchunk} to the desired (integer) chunk size; to disable 
  chunking set \code{nchunk=NA}.
  
  Both \code{MAP} and \code{REDUCE} are applied in the distributed 
  step (``on the worker``). There is no built-in ability to combine 
  results across workers in the distributed step.
}

\value{
  \itemize{
    \item{reduceByFile:}{
      When \code{summarize=FALSE} the return value is a \code{list} or
      the value from the final invocation of \code{REDUCE}. When
      \code{summarize=TRUE} output is a \code{SummarizedExperiment}.
      When \code{ranges} is a \code{GenomicFiles} object data from 
      \code{rowData}, \code{colData} and \code{exptData} are transferred 
      to the \code{SummarizedExperiment}.
    }
    \item{reduceFiles:}{
      A \code{list} or the value returned by the final invocation of
      \code{REDUCE}.
    }
  }
}

\seealso{
  \itemize{
    \item \link{reduceRanges}
    \item \link{reduceByRange}
    \item \link{GenomicFiles-class}
  }
}

\author{
  Martin Morgan <mtmorgan@fhcrc.org> and 
  Valerie Obenchain <vobencha@fhcrc.org>
}

\examples{
## This example counts junction reads for ranges of interest in BAM files.

if (all(require(RNAseqData.HNRNPC.bam.chr14) &&
        require(GenomicAlignments))) {
  fls <- RNAseqData.HNRNPC.bam.chr14_BAMFILES  ## 8 bam files
  
  ## Ranges of interest.
  gr <- GRanges("chr14", IRanges(c(19100000, 106000000), width=1e7))
  
  ## MAP outputs a table of junction counts per range.
  MAP <- function(range, file, ...) {
      library(GenomicAlignments)
      param = ScanBamParam(which=range)
      gal = readGAlignments(file, param=param)
      table(njunc(gal))
  } 

  ## -----------------------------------------------------------------------
  ## reduceByFile:
  ## -----------------------------------------------------------------------

  ## reduceByFile() with a MAP and no REDUCE returns junction counts for
  ## each range / file combination.
  counts1 <- reduceByFile(gr, fls, MAP)
  length(counts1)          ## 8 files
  elementLengths(counts1)  ## 2 ranges each
  
  ## Each element contains a table of counts for each range.
  counts1[[1]]

  ## This reducer sums the number of records in each range 
  ## with exactly 1 junction.
  REDUCE <- function(mapped, ...)
      sum(sapply(mapped, "[", "1"))
  
  reduceByFile(gr, fls, MAP, REDUCE)

  ## -----------------------------------------------------------------------
  ## reduceFiles:
  ## -----------------------------------------------------------------------

  ## This function treats all ranges as a group and is useful when 
  ## differentiating between ranges is not important. In the 
  ## reduceByFiles() example junctions were counted for each range 
  ## individually which allowed us to see results for the individual
  ## ranges and combine them on the fly in a specific way. reduceFiles() 
  ## counts junctions for all ranges simultaneously. Any further manipulations
  ## must be done to the result as a whole.
  counts2 <- reduceFiles(gr, fls, MAP)

  ## reduceFiles() returns counts for all ranges.
  counts2[[1]]
  ## reduceByFile() returns counts for each range separately.
  counts1[[1]]

  ## -----------------------------------------------------------------------
  ## Methods for the GenomicFiles class:
  ## -----------------------------------------------------------------------
 
  ## Both reduceByFiles() and reduceFiles() can operate on a GenomicFiles
  ## object. The class is a light-weight form of SummarizedExperiment
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
  
  ## Results are in the assays slot named 'data'.
  assays(se) 
}
}

\keyword{methods}
