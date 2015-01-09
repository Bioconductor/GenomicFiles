\name{reduceByFile}

\alias{reduceByFile}
\alias{reduceByFile,GRanges,ANY-method}
\alias{reduceByFile,GRangesList,ANY-method}
\alias{reduceByFile,GenomicFiles,missing-method}

\alias{reduceFiles}

\title{Parallel computations by files}

\description{
  Computations are distributed in parallel by file. Data subsets are 
  extracted and manipulated (MAP) and optionally combined (REDUCE)
  within a single file.
} 

\usage{
\S4method{reduceByFile}{GRanges,ANY}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)
\S4method{reduceByFile}{GRangesList,ANY}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)
\S4method{reduceByFile}{GenomicFiles,missing}(ranges, files, MAP, 
    REDUCE, ..., summarize=FALSE, iterate=TRUE, init)

reduceFiles(ranges, files, MAP, REDUCE, ..., init)
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
    is set to FALSE. This argument applies to \code{reduceByFile} only
    (\code{reduceFiles} calls MAP a single time on each worker).

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

if (all(require(RNAseqData.HNRNPC.bam.chr14) &&
        require(GenomicAlignments))) {

  ## -----------------------------------------------------------------------
  ## Count junction reads in BAM files
  ## -----------------------------------------------------------------------
  fls <- RNAseqData.HNRNPC.bam.chr14_BAMFILES  ## 8 bam files
 
  ## Ranges of interest.
  gr <- GRanges("chr14", IRanges(c(19100000, 106000000), width=1e7))
 
  ## MAP outputs a table of junction counts per range.
  MAP <- function(range, file, ...) {
      library(GenomicAlignments)   ## for readGAlignments(), ScanBamParam()
      param = ScanBamParam(which=range)
      gal = readGAlignments(file, param=param)
      table(njunc(gal))
  } 

  ## -----------------------------------------------------------------------
  ## reduceByFile:

  ## With no REDUCE, counts are computed for each range / file combination.
  counts1 <- reduceByFile(gr, fls, MAP)
  length(counts1)          ## 8 files
  elementLengths(counts1)  ## 2 ranges each
 
  ## Tables of counts for each range:
  counts1[[1]]

  ## With a REDUCE, results are combined on the fly. This reducer sums the 
  ## number of records in each range with exactly 1 junction.
  REDUCE <- function(mapped, ...)
      sum(sapply(mapped, "[", "1"))
 
  reduceByFile(gr, fls, MAP, REDUCE)

  ## -----------------------------------------------------------------------
  ## reduceFiles:

  ## All ranges are treated as a single group:
  counts2 <- reduceFiles(gr, fls, MAP)

  ## Counts are for all ranges grouped:
  counts2[[1]]

  ## Contrast the above with that from reduceByFile() where counts 
  ## are for each range separately:
  counts1[[1]]

  ## -----------------------------------------------------------------------
  ## Methods for the GenomicFiles class:
 
  ## Both reduceByFiles() and reduceFiles() can operate on a GenomicFiles
  ## object.
  colData <- DataFrame(method=rep("RNASeq", length(fls)),
                       format=rep("bam", length(fls)))
  gf <- GenomicFiles(files=fls, rowData=gr, colData=colData)
  gf
  
  ## Subset on ranges or files for different experimental runs.
  dim(gf)
  gf_sub <- gf[2, 3:4]
  dim(gf_sub)
  
  ## When summarize = TRUE and no REDUCE is given, the output is a 
  ## SummarizedExperiment object.
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
