\name{FaFileViews}
\docType{class}

% Class:
\alias{FaFileViews-class}
\alias{FaFileViews}
\alias{FaFileViews,character-method}
\alias{FaFileViews,FaFileList-method}
\alias{FaFileViews,missing-method}

\title{Views into a set of Fasta files}

\description{
  This class is deprecated. Use \code{GenomicFiles} instead.
  Use \code{FaFileViews()} to reference a set of disk-based Fasta files to be
  processed
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{FaFileViews(fileList,
        fileSample=DataFrame(row.names=make.unique(basename(path(fileList)))),
        fileRange, fileExperiment = list(), yieldSize = NA_integer_, 
        .views_on_file = new.env(parent=emptyenv()), ...)}:
      This constructor is a generic function with dispatch on argument
      \code{fileList}. Methods exist for \link{FaFileList} and character 
      (vector of file names). 
    }
  }
}

\section{Accessors}{
  All accessor-like methods defined for \code{GenomicFileViews} objects
  work on \code{FaFileViews} objects. See ?\code{GenomicFileViews} for details.
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
      Parallel computations are distributed by files in \code{fileList}
      with the option to provide MAP and REDUCE functions across ranges 
      and / or files.
    }
    \item{reduceByRange}{
      Parallel computations are distributed by ranges in \code{fileRange}
      with the option to provide MAP and REDUCE functions across ranges 
      and / or files.
    }
  }
}

\section{Arguments}{
  \describe{
    \item{}{
       \code{fileList}:
         A character() vector of Fasta path names or a 
         \code{\linkS4class{FaFileList}}.
    }
    \item{}{
      \code{fileSample}:
        A \link[S4Vectors]{DataFrame} instance with as many rows as 
        \code{length(fileList)}, containing sample information associated 
        with each path.
    }
    \item{}{
      \code{fileRange}:
        A \code{\link[GenomicRanges]{GRanges}}, or missing instance 
        with ranges defined on the spaces of the Fasta files. Ranges are 
        \emph{not} validated against  the Fasta files.
    }
    \item{}{
      \code{fileExperiment}:
        A list() containing additional information about the experiment.
    }
    \item{}{
      \code{yieldSize}: An integer specifying number of records to process
    }
    \item{}{
        \code{.views_on_file}: An environment; currently under development
    }
    \item{}{
      \code{...}: Additional arguments.
    }
    \item{}{
      \code{x}, \code{object}: An instance of \code{FaFileViews}.
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
      \code{file}: An instance of \code{FaFileViews}.
    }
    \item{}{
      \code{index}: Not used.
    }
    \item{}{
      \code{param}:
      Unused option for \code{FaFileViews} object. \code{fileRange}
      are used to specify ranges to query.
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
    \item \link{GenomicFileViews-class}.
  }
}

\author{
  Martin Morgan <mtmorgan@fhcrc.org.> and 
  Valerie Obenchain <vobencha@fhcrc.org>
}

\examples{
fa <- system.file("extdata", "ce2dict1.fa", package="Rsamtools",
                  mustWork=TRUE)
ffv <- FaFileViews(c(fa, fa), fileRange=scanFaIndex(fa))
ffv
}

\keyword{classes}
\keyword{methods}
