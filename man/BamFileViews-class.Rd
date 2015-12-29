\name{BamFileViews}
\docType{class}

% Class:
\alias{BamFileViews-class}
\alias{BamFileViews}
\alias{BamFileViews,character-method}
\alias{BamFileViews,BamFileList-method}
\alias{BamFileViews,missing-method}
\alias{BamFileViews,ANY-method}

% Methods:
\alias{countBam,BamFileViews-method}
\alias{scanBam,BamFileViews-method}
\alias{summarizeOverlaps,GRanges,BamFileViews-method}
\alias{summarizeOverlaps,GRangesList,BamFileViews-method}

\title{Views into a set of BAM files}

\description{
  This class is defunct. Use \code{GenomicFiles} instead.
  Use \code{BamFileViews()} to reference a set of disk-based BAM files to be
  processed 
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{BamFileViews(fileList,
        fileSample=DataFrame(row.names=make.unique(basename(path(fileList)))),
        fileRange, fileExperiment = list(), yieldSize = NA_integer_, 
        .views_on_file = new.env(parent=emptyenv()), ...)}:
      This constructor is a generic function with dispatch on argument
      \code{fileList}. Methods exist for \link{BamFileList} and character 
      (vector of file names). 
    }
  }
}

\section{Accessors}{
  All accessor-like methods defined for \code{GenomicFileViews} objects
  work on \code{BamFileViews} objects. See ?\code{GenomicFileViews} for details.
  \itemize{
    \item{fileList(x); fileList(x) <- value}
    \item{fileSample(x); fileSample(x) <- value}
    \item{fileRange(x); fileRange(x) <- value}
    \item{fileExperiment(x); fileExpermient(x) <- value}
    \item{yieldSize(x); yieldSize(x) <- value}
  }
}

\section{Methods}{
  \describe{
    \item{}{
      \code{"["}: Subset the object by \code{fileRange} or \code{fileSample}.
    }
    \item{reduceByFile}{
      Computations are distributed in parallel by files in \code{fileList}
      with the option to provide MAP and REDUCE functions across ranges 
      and / or files.
    }
    \item{reduceByRange}{
      Computations are distributed in parallel by ranges in \code{fileRange}
      with the option to provide MAP and REDUCE functions across ranges 
      and / or files.
    }
    \item{summarizeOverlaps}{
      Computations are distributed in parallel by files in \code{fileList}.
      Ranges in the \code{fileRange} slot take precedence over ranges
      in \code{param}. The return value is a \code{SummarizedExperiment} 
      object.
    }
    \item{countBam}{
      Computations are distributed in parallel by files in \code{fileList}.
      Ranges in the \code{fileRange} slot take precedence over ranges
      in \code{param}. The return value is a list of \code{data.frame}s,
      one per file.
    }
    \item{scanBam}{
      Computations are distributed in parallel by files in \code{fileList}.
      Ranges in the \code{fileRange} slot take precedence over ranges
      in \code{param}. The return value is a list of \code{list}s, one per file.
    }
  }
}

\section{Arguments}{
  \describe{
    \item{}{
       \code{fileList}:
         A character() vector of BAM path names or a 
         \code{\linkS4class{BamFileList}}.
    }
    \item{}{
      \code{fileSample}:
        A \link[S4Vectors]{DataFrame} instance with as many rows as 
        \code{length(fileList)}, containing sample information associated 
        with each path.
    }
    \item{}{
      \code{fileRange}:
        A \code{\link[GenomicRanges]{GRanges}}, or missing instance with 
        ranges defined on the spaces of the BAM files. Ranges are \emph{not} 
        validated against the BAM files.
    }
    \item{}{
      \code{fileExperiment}:
        A list() containing additional information about the experiment.
    }
    \item{}{
      \code{yieldSize}: An integer specifying number of records to process
    }
    \item{}{
        \code{.views_on_file}: An enviornment; currently under development
    }
    \item{}{
      \code{...}: Additional arguments.
    }
    \item{}{
      \code{x}, \code{object}: An instance of \code{BamFileViews}.
    }
    \item{}{
      \code{value}: An object of appropriate type to replace content.
    }
    \item{}{
      \code{i}: 
      During subsetting, a logical or numeric index into \code{fileRange}.
    }
    \item{}{
      \code{j}: 
      During subsetting, a logical or numeric index into \code{fileSample} 
      and \code{fileList}.
    }
    \item{}{
      \code{file}: An instance of \code{BamFileViews}.
    }
    \item{}{
      \code{index}: Not used.
    }
    \item{}{
      \code{param}:
      An optional \code{\linkS4class{ScanBamParam}} instance to further 
      influence scanning or counting.
    }
  }
}

\section{Slots}{
  \describe{
  Inherited from \code{GenomicFileViews} class:
    \itemize{
      \item fileList
      \item fileSample
      \item fileRange
      \item fileExperiment
      \item yieldSize
      \item .views_on_file
    }
  }
}

\seealso{
  \itemize{
    \item \link{GenomicFileViews-class} class.
    \item \link{reduceByFile} and \link{reduceByRange} methods.
  }
}

\author{
  Martin Morgan and Valerie Obenchain
}

\examples{
## See ?GenomicFiles.
}

\keyword{classes}
\keyword{methods}
