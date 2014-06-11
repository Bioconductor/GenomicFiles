\name{GenomicFiles}
\docType{class}

% Class:
\alias{GenomicFiles}
\alias{class:GenomicFiles}
\alias{GenomicFiles-class}

% Constructors:
\alias{GenomicFiles,character-method}
\alias{GenomicFiles,list-method}
\alias{GenomicFiles,List-method}
\alias{GenomicFiles,missing-method}

% Accessors:
\alias{files}
\alias{files<-}
\alias{files,GenomicFiles-method}
\alias{files<-,GenomicFiles,ANY-method}
\alias{colData,GenomicFiles-method}
\alias{colData<-,GenomicFiles,DataFrame-method}
\alias{rowData,GenomicFiles-method}
\alias{rowData<-,GenomicFiles,DataFrame-method}
\alias{exptData,GenomicFiles-method}
\alias{exptData<-,GenomicFiles,DataFrame-method}

% Methods:
\alias{dim,GenomicFiles-method}
\alias{dimnames,GenomicFiles-method}
\alias{dimnames<-,GenomicFiles,ANY-method}
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

  Methods set on \code{GenomicFiles} return a \code{SummarizedExperiment}
  where data from \code{rowData}, \code{colData} and \code{exptData}
  are transferred from one class to the other.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{GenomicFiles(files, rowData=DataFrame(), colData=DataFrame(),
                         exptData=SimpleList(), ...)}:
    }
  }
}

\section{Accessors}{
  In the code below, \code{x} is a GenomicFiles object.

  \describe{
    item{files(x), files(x) <- value}{
      Get or set the files on \code{x}. \code{value} can be 
      a character() of file paths or a List of file objects
      such as BamFile, BigWigFile etc.
    }
    \item{rowData, rowData(x) <- value}{
      Get or set the rowData on \code{x}. \code{value} must be
      a \code{DataFrame} instance with ranges or indices defined 
      on the spaces (genomic position) of the files.
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
%    \item{reduceByFile}{
%      Parallel computations are distributed by files in \code{fileList}
%      with the option to provide MAP and REDUCE functions across ranges 
%      and / or files.
%    }
%    \item{reduceByRange}{
%      Parallel computations are distributed by ranges in \code{fileRange}
%      with the option to provide MAP and REDUCE functions across ranges 
%      and / or files.
%    }
  }
}

\seealso{
  \itemize{
    \item \link{reduceByFile} and \link{reduceByRange} methods.
  }
}

\author{
  Valerie Obenchain <vobencha@fhcrc.org>
}

\examples{
library(RNAseqData.HNRNPC.bam.chr14)
fls <- RNAseqData.HNRNPC.bam.chr14_BAMFILES
colData <- DataFrame(method=rep("RNASeq", length(fls)), 
                     format=rep("bam", length(fls)))
rowData <- DataFrame(GRanges("chr14", IRanges(2e5*1:3, width=1e5)))
gf <- GenomicFiles(fls, rowData=rowData, colData=colData)
gf 
}

\keyword{classes}
\keyword{methods}
