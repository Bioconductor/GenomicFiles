\name{GenomicFiles}
\docType{class}

% Class:
\alias{GenomicFiles}
\alias{class:GenomicFiles}
\alias{GenomicFiles-class}

% Constructors:
\alias{GenomicFiles,GenomicRanges_OR_GRangesList,character-method}
\alias{GenomicFiles,GenomicRanges_OR_GRangesList,List-method}
\alias{GenomicFiles,GenomicRanges_OR_GRangesList,list-method}
\alias{GenomicFiles,missing,ANY-method}
\alias{GenomicFiles,missing,missing-method}

% Accessors:
\alias{files<-}
\alias{files}
\alias{files,GenomicFiles-method}
\alias{files<-,GenomicFiles,character-method}
\alias{files<-,GenomicFiles,List-method}
\alias{dimnames<-,GenomicFiles,list-method}
\alias{colData<-,GenomicFiles,DataFrame-method}

% Methods:
\alias{[,GenomicFiles,ANY,ANY-method}
\alias{[,GenomicFiles,ANY,ANY,ANY-method}
\alias{show,GenomicFiles-method}

\title{GenomicFiles objects}

\description{
  The \code{GenomicFiles} class is a matrix-like container where rows
  represent ranges of interest and columns represent files. The class 
  is designed for byFile or byRange queries.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicFiles(rowRanges, files, colData=DataFrame(),
                         metadata=list(), ...)}:
    }
  }
}

\section{Details}{
    \code{GenomicFiles} inherits from the \code{RangedSummarizedExperiment}
    class in the \code{SummarizedExperiment} package. Currently, no use
    is made of the \code{elementMetadat} and \code{assays} slots. This may
    change in the future.
}

\section{Accessors}{
  In the code below, \code{x} is a GenomicFiles object.

  \describe{
    \item{rowRanges, rowRanges(x) <- value}{
      Get or set the rowRanges on \code{x}. \code{value} can be
      a \code{GRanges} or \code{GRangesList} representing ranges or indices 
      defined on the spaces (position) of the files.
    }
    \item{files(x), files(x) <- value}{
      Get or set the files on \code{x}. \code{value} can be 
      a character() of file paths or a List of file objects
      such as BamFile, BigWigFile, FaFile, etc.
    }
    \item{colData, colData(x) <- value}{
      Get or set the colData on \code{x}. \code{value} must be
      a \code{DataFrame} instance describing the files. The number 
      of rows must match the number of files. Row names, if present, 
      become the column names of the \code{GenomicFiles}.
    }
    \item{metadata, metadata(x) <- value}{
      Get or set the metadata on \code{x}. \code{value} must be
      a SimpleList of arbitrary content describing the overall 
      experiment.
    }
    \item{dimnames, dimnames(x) <- value}{
      Get or set the row and column names on \code{x}.
    }
  }
}

\section{Methods}{
  In the code below, \code{x} is a GenomicFiles object.

  \describe{
    \item{[}{
      Subset the object by \code{fileRange} or \code{fileSample}.
    }
    \item{show}{
      Compactly display the object.
    }
    \item{reduceByFile}{
     Extract, manipulate and combine data defined in \code{rowRanges} 
     within the files specified in \code{files}. See ?\code{reduceByFile} 
     for details.
    }
    \item{reduceByRange}{
     Extract, manipulate and combine data defined in \code{rowRanges} 
     across the files specified in \code{files}. See ?\code{reduceByRange} 
     for details.
    }
  }
}

\seealso{
  \itemize{
    \item \link{reduceByFile} and \link{reduceByRange} methods.
    \item \link[SummarizedExperiment]{SummarizedExperiment} objects in the 
          \pkg{SummarizedExperiment} package.
  }
}

\author{
  Martin Morgan and Valerie Obenchain
}

\examples{
## -----------------------------------------------------------------------
## Basic Use
## -----------------------------------------------------------------------

if (require(RNAseqData.HNRNPC.bam.chr14)) { 
  fl <- RNAseqData.HNRNPC.bam.chr14_BAMFILES
  rd <- GRanges("chr14", 
                 IRanges(c(62262735, 63121531, 63980327), width=214700))
  cd <- DataFrame(method=rep("RNASeq", length(fl)),
                  format=rep("bam", length(fl)))

  ## Construct an instance of the class:
  gf <- GenomicFiles(files = fl, rowRanges = rd, colData = cd) 
  gf

  ## Subset on ranges or files for different experimental runs.
  dim(gf)
  gf_sub <- gf[2, 3:4]
  dim(gf_sub)
 
  ## When summarize = TRUE and no REDUCE is provided the reduceBy* 
  ## functions output a SummarizedExperiment object.
  MAP <- function(range, file, ...) {
      requireNamespace("GenomicFiles", quietly=TRUE) ## for coverage()
      requireNamespace("Rsamtools", quietly=TRUE)     ## for ScanBamParam()
      param = Rsamtools::ScanBamParam(which=range)
      GenomicFiles::coverage(file, param=param)[range]
  } 
  se <- reduceByRange(gf, MAP=MAP, summarize=TRUE)
  se
 
  ## Data from the rowRanges, colData and metadata slots in the
  ## GenomicFiles are transferred to the SummarizedExperiment.
  colData(se)
 
  ## Results are in the assays slot.
  assays(se) 
}
 
## -----------------------------------------------------------------------
## Managing cached or remote files with GenomicFiles
## -----------------------------------------------------------------------

## The GenomicFiles class can manage cached or remote files and their 
## associated ranges.

\dontrun{
## Files from AnnotationHub can be downloaded and cached locally.
library(AnnotationHub)
hub = AnnotationHub()
hublet = query(hub, c("files I'm", "interested in"))
# cache (if need) and return local path to files
fls = cache(hublet)

## An alternative to the local file paths is to use urls to a remote file.
## This approach could be used with something like rtracklayer::bigWig which
## supports remote file queries.
urls = hublet$sourceurls

## Define ranges of interest and use GenomicFiles to manage.
rngs = GRanges("chr10", IRanges(c(100000, 200000), width=1))
gf = GenomicFiles(rngs, fls)

## As an example, one could create a matrix from data extracted
## across multiple BED files.
MAP = function(rng, fl) {
    requireNamespace("rtracklayer", quietly=TRUE)  ## import, BEDFile
    rtracklayer::import(rtracklayer::BEDFile(fl), which=rng)$name
}
REDUCE = unlist
xx = reduceFiles(gf, MAP=MAP, REDUCE=REDUCE)
mcols(rngs) = simplify2array(xx)

## Data and ranges can be stored in a SummarizedExperiment.
SummarizedExperiment(list(my=simplify2array(xx)), rowRanges=rngs)
}
}

\keyword{classes}
\keyword{methods}
