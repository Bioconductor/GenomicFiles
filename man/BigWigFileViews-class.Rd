\name{BigWigFileViews}
\docType{class}

% Class:
\alias{BigWigFileViews-class}
\alias{BigWigFileViews}
\alias{BigWigFileViews,BigWigFileList-method}
\alias{BigWigFileViews,character-method}
\alias{BigWigFileViews,missing-method}

% Methods:
\alias{coverage,BigWigFileViews-method}
\alias{summary,BigWigFileViews-method}

\title{Views into a set of BigWig files}

\description{
  This class is deprecated. Use \code{GenomicFiles} instead.
  Use \code{BigWigFileViews()} to reference a set of disk-based BigWig files 
  to be processed (e.g., queried using \code{\link[IRanges]{coverage}}) as 
  a single \sQuote{experiment}.
}

\section{Constructor}{
  \describe{
    \item{}{
      \code{BigWigFileViews(fileList,
        fileSample=DataFrame(row.names=make.unique(basename(path(fileList)))),
        fileRange, fileExperiment = list(), yieldSize = NA_integer_, 
        .views_on_file = new.env(parent=emptyenv()), ...)}:
      This constructor is a generic function with dispatch on argument
      \code{fileList}. Methods exist for 
      \code{\link[rtracklayer]{BigWigFileList}} and \code{character()} 
      (vector of file names). 
    }
  }
}

\section{Accessors}{
  All accessor-like methods defined for \code{GenomicFileViews} objects
  work on \code{BigWigFileViews} objects. See ?\code{GenomicFileViews} for details.
  \itemize{
    \item{fileList(x); fileList(x) <- value}
    \item{fileSample(x); fileSample(x) <- value}
    \item{fileRange(x); fileRange(x) <- value}
    \item{fileExperiment(x); fileExpermient(x) <- value}
    \item{yieldSize(x); yieldSize(x) <- value}
  }
}

\section{Subsetting}{
  \describe{
    \item{}{
      \code{"["}: Subset the object by \code{fileRange} or \code{fileSample}.
    }
  }
}

\section{Other methods}{
  In the code snippets below, \code{x} and \code{object} are 
  BigWigFileViews objects.

  \describe{
    \item{}{
       \code{coverage(x, ..., by = "file", summarize = TRUE,
                      as = "RleList")}:
       Computes coverage with the \code{import} function from
       \code{rtracklayer} for each file in \code{fileList(x)} and 
       each range in \code{fileRange(x)}. Work is divided
       in parallel as specified in the \code{by} argument. Results
       are returned as a list unless \code{summarize=TRUE} in
       which case the data are in the \code{assays} slot of
       a SummarizedExperiment object. Data type are controlled with
       the \code{as} argument. See ?\code{'import,BigWigFile-method'}
       for details.
    }
    \item{}{
       \code{summary(object, ..., by = "file", summarize = TRUE)}:
       Computes summary statistics with the \code{summary} function
       from \code{rtracklayer} for each file in \code{fileList(object)} 
       and each range in \code{fileRange(object)}. Work is divided in 
       parallel as specified in the \code{by} argument. Results are 
       returned as a list unless \code{summarize=TRUE} in which case the 
       data are in the \code{assays} slot of a SummarizedExperiment 
       object. Summary statistics are controlled with the \code{type}
       argument passed to \code{summary}. See
       ?\code{'summary,BigWigFile-method'} for details.
    }
  }
}

\section{Arguments}{
  \describe{
    \item{}{
       \code{fileList}:
         A character() vector of BigWig path names or a 
         \code{\link[rtracklayer]{BigWigFileList}}.
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
        ranges defined on the spaces of the BigWig files. Ranges are 
        \emph{not} validated against the BigWig files.
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
      \code{x}, \code{object}: An instance of \code{BigWigFileViews}.
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
      \code{file}: An instance of \code{BigWigFileViews}.
    }
    \item{}{
      \code{index}: Not used.
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
  Michael Love <michaelisaiahlove@gmail.com>,
  Valerie Obenchain <vobencha@fhcrc.org>,
  Martin Morgan <mtmorgan@fhcrc.org.>
}


\examples{
if (.Platform$OS.type != "windows") {

  register(SerialParam())

  ## ---------------------------------------------------------------------
  ## BigWigFileView Objects
  ## ---------------------------------------------------------------------
  fl <- system.file("tests", "test.bw", package = "rtracklayer")
  gr <- GRanges(Rle(c("chr2", "chr19"), c(4, 2)),
                IRanges(1 + c(200, 250, 500, 550, 1450, 1750), width=100))
  bwfv <- BigWigFileViews(c(fl, fl), fileRange=gr)
 
  ## Object dimensions are range by sample (files).
  bwfv 
  fileList(bwfv)
  fileRange(bwfv)
 
  ## Subset by range(s) or sample(s).
  bwfv[1:3,]
  bwfv[1:3,2]

  ## ---------------------------------------------------------------------
  ## coverage()
  ## ---------------------------------------------------------------------
  ## coverage() on a 'BigWigFileViews' object returns coverage for all 
  ## files in fileList() for the ranges in fileRange().
 
  ## When 'summarize = FALSE' results are returned as a list. The
  ## 'by' argument specifies grouping by range or by file.
  cv1 <- coverage(bwfv, by = "file", summarize = FALSE)
  length(cv1)
  elementLengths(cv1)
 
  cv2 <- coverage(bwfv, by = "range", summarize = FALSE)
  length(cv2)
  elementLengths(cv2)
 
  ## When 'summarize = TRUE' output is a 'SummarizedExperiment'.
  cv3 <- coverage(bwfv, summarize = TRUE)
  assays(cv3)
 
  ## FIXME: improve show method for 'Lists' in a matrix
  ## FIXME: assays as 'List' should respect 'withDimnames'
  assays(cv3, withDimnames=TRUE)[[1]]
 
  ## ---------------------------------------------------------------------
  ## summary() 
  ## ---------------------------------------------------------------------
  ## summary() on a 'BigWigFileViews' object returns the 'type' statistic
  ## for all files in fileList() for the ranges in fileRange().

  sm1 <- summary(bwfv, type = "mean", by = "file", summarize = TRUE) 
  sm2 <- summary(bwfv, type = "max", by = "range", summarize = FALSE) 
} 
}

\keyword{classes}
\keyword{methods}
