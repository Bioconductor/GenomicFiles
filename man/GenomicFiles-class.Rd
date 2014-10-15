\name{GenomicFiles}
\docType{class}

% Class:
\alias{GenomicFiles}
\alias{class:GenomicFiles}
\alias{GenomicFiles-class}

% Constructors:
\alias{GenomicFiles,GenomicRangesORGRangesList,character-method}
\alias{GenomicFiles,GenomicRangesORGRangesList,List-method}
\alias{GenomicFiles,GenomicRangesORGRangesList,list-method}
\alias{GenomicFiles,missing,ANY-method}
\alias{GenomicFiles,missing,missing-method}

% Accessors:
\alias{files}
\alias{files<-}
\alias{files,GenomicFiles-method}
\alias{files<-,GenomicFiles,character-method}
\alias{files<-,GenomicFiles,List-method}
\alias{colData,GenomicFiles-method}
\alias{colData<-,GenomicFiles,DataFrame-method}
\alias{rowData,GenomicFiles-method}
\alias{rowData<-,GenomicFiles,GRangesORGRangesList-method}
\alias{exptData,GenomicFiles-method}
\alias{exptData<-,GenomicFiles,DataFrame-method}

% Methods:
\alias{dim,GenomicFiles-method}
\alias{dimnames,GenomicFiles-method}
\alias{dimnames<-,GenomicFiles,NULL-method}
\alias{dimnames<-,GenomicFiles,list-method}
\alias{[,GenomicFiles,ANY,missing-method}
\alias{[,GenomicFiles,missing,ANY-method}
\alias{[,GenomicFiles,ANY,ANY-method}
\alias{[,GenomicFiles,ANY,missing-method}
\alias{show,GenomicFiles-method}

\title{GenomicFiles objects}

\description{
  The \code{GenomicFiles} class is a matrix-like container where rows
  represent ranges of interest and columns represent files. The class 
  is designed for byFile or byRange queries and can be thought of as 
  a precursor (or skeleton) of the \code{SummarizedExperiment} 
  class.

  \code{GenomicFiles} has slots for \code{files}, \code{rowData}, 
  \code{colData} and \code{exptData}. In contrast to 
  \code{SummarizedExperiment}, \code{GenomicFiles} does not contain an 
  \code{assays} slot and the \code{rowData} is a \code{DataFrame} 
  instead of a \code{GenomicRanges}. 
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicFiles(rowData, files, colData=DataFrame(),
                         exptData=SimpleList(), ...)}:
    }
  }
}

\section{Accessors}{
  In the code below, \code{x} is a GenomicFiles object.

  \describe{
    \item{rowData, rowData(x) <- value}{
      Get or set the rowData on \code{x}. \code{value} can be
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
    \item{exptData, exptData(x) <- value}{
      Get or set the exptData on \code{x}. \code{value} must be
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
     Extract, manipulate and combine data defined in \code{rowData} 
     within the files specified in \code{files}. See ?\code{reduceByFile} 
     for details.
    }
    \item{reduceByRange}{
     Extract, manipulate and combine data defined in \code{rowData} 
     across the files specified in \code{files}. See ?\code{reduceByRange} 
     for details.
    }
  }
}

\seealso{
  \itemize{
    \item \link{reduceByFile} and \link{reduceByRange} methods.
  }
}

\author{
  Martin Morgan <mtmorgan@fhcrc.org> and 
  Valerie Obenchain <vobencha@fhcrc.org>
}

\examples{
if (require(RNAseqData.HNRNPC.bam.chr14)) { 
  fls <- RNAseqData.HNRNPC.bam.chr14_BAMFILES
  gr <- GRanges("chr14", IRanges(c(62262735, 63121531, 63980327),
                width=214700))
  
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
  MAP <- function(range, file, ...) {
      library(GenomicRanges)
      param = ScanBamParam(which=range)
      coverage(file, param=param)[range]
  } 
  se <- reduceByRange(gf, MAP=MAP, summarize=TRUE)
  se
  
  ## Data from the rowData, colData and exptData slots in the
  ## GenomicFiles are transferred to the SummarizedExperiment.
  colData(se)
  
  ## Results are in the assays slot.
  assays(se) 
}
}

\keyword{classes}
\keyword{methods}
