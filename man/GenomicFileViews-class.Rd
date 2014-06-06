\name{GenomicFileViews}
\docType{class}

% Class:
\alias{GenomicFileViews}
\alias{class:GenomicFileViews}
\alias{GenomicFileViews-class}

% Accessors:
\alias{fileList}
\alias{fileList<-}
\alias{fileSample}
\alias{fileSample<-}
\alias{fileRange}
\alias{fileRange<-}
\alias{fileExperiment}
\alias{fileExperiment<-}

% Methods:
\alias{names,GenomicFileViews-method}
\alias{names<-,GenomicFileViews-method}
\alias{dim,GenomicFileViews-method}
\alias{dimnames,GenomicFileViews-method}
\alias{dimnames<-,GenomicFileViews,ANY-method}
\alias{[,GenomicFileViews,ANY,missing-method}
\alias{[,GenomicFileViews,missing,ANY-method}
\alias{[,GenomicFileViews,ANY,ANY-method}
\alias{[,GenomicFileViews,ANY,missing-method}
\alias{yieldSize,GenomicFileViews-method}
\alias{yieldSize<-,GenomicFileViews-method}
\alias{show,GenomicFileViews-method}

\title{Views into a set of files}

\description{
  \code{GenomicFileViews} is a VIRTUAL class used to reference a set of 
  disk-based files to be queried across views (ranges).
}

\section{Objects from the Class}{
  \code{GenomicFileViews} is a VIRTUAL class not intended for instantiation by
  the user. The class serves as a parent for concrete subclasses such 
  as \code{BamFileViews}, \code{FaFileViews}, \code{TabixFileViews} etc.
}

\section{Slots}{

  \describe{
    \item{fileList}{
      List of of length >= 2 containing the file path and index names. 
      List names must include \sQuote{path} and \sQuote{index}.
    }
    \item{fileSample}{
        A \code{\link[IRanges]{DataFrame}} instance with as many rows as 
        \code{length(fileList)}, containing sample information associated 
        with each path.
    }
    \item{fileRange}{
      A \code{\link[GenomicRanges]{GRanges}} instance with ranges defined 
      on the spaces (genomic position) of the files.
    }
    \item{fileExperiment}{
      A \code{list} containing additional information about the experiment.
    }
    \item{yieldSize}{
        An integer specifying the data chunk size.
    }
    \item{.views_on_file}{
        An environment. Under construction / future use.
    }
  }
}

\section{Accessors}{
  In the code snippets below, \code{x} is a GenomicFileViews object.

  \describe{
    item{fileList(x), fileList(x) <- value}{
      Get or set the fileList on \code{x}. \code{value} must be 
      a List with list elements appropriate for the subclass.
    }
    \item{fileSample, fileSample(x) <- value}{
      Get or set the fileSample on \code{x}. \code{value} must be
      a \code{\link[IRanges]{DataFrame}} instance with as many rows as 
      \code{length(fileList)}, containing sample information associated 
      with each file.
    }
    \item{fileRange, fileRange(x) <- value}{
      Get or set the fileSample on \code{x}. \code{value} must be
      a \code{\link[GenomicRanges]{GRanges}} instance.
    }
    \item{fileExperiment, fileExperiment(x) <- value}{
      Get or set the fileExperiment on \code{x}. \code{value} must be
      a list().
    }
    \item{yieldSize, yieldSize(x) <- value}{
      Get or set the yieldSize on \code{x}. \code{value} must be
      an integer.
    }
    \item{names, names(x) <- value}{
      Get or set the names on \code{x}. These are the column
      names of the \code{GenomicFileViews} instance corresponding to
      the paths in fileList.
    }
    \item{dimnames, dimnames(x) <- value}{
      Get or set the row and column names on \code{x}.
    }
  }
}

\section{Methods}{
  In the code snippets below, \code{x} is a GenomicFileViews object.

  \describe{
    \item{[}{
      Subset the object by \code{fileRange} or \code{fileSample}.
    }
    \item{show}{
      Compactly display the object.
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

\seealso{
  \itemize{
    \item \link{BamFileViews-class} class.
    \item \link{BigWigFileViews-class} class.
    \item \link{reduceByFile} and \link{reduceByRange} methods.
  }
}

\author{
  Martin Morgan <mtmorgan@fhcrc.org.> and 
  Valerie Obenchain <vobencha@fhcrc.org>
}

\examples{
fl <- system.file("extdata", "ex1.bam", package="Rsamtools",
                  mustWork=TRUE) 
bfv <- BamFileViews(fl, fileRange = GRanges("seq1", IRanges(1, 10)))
showClass(class(bfv))
 
## See ?BamFileViews, ?BigWigFileViews, ?reduceBy for additional examples.
}

\keyword{classes}
\keyword{methods}
