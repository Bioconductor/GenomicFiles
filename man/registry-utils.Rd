\name{registry-utils}

\alias{registry-utils}

\alias{registerFileType}
\alias{findTypeRegistry}
\alias{makeFileType}

\title{Functions for creating and searching a registry of file types.}

\description{
  Functions for creating and searching a registry of file types
  based on file extension.
}

\usage{
  registerFileType(type, package, regex)
  findTypeRegistry(fnames)
  makeFileType(fnames, ..., regex=findTypeRegistry(fnames))
}

\arguments{
  \item{type}{
    The List class the file is associated with such as
    BamFileList, BigWigFileList, FaFileList.
  }
  \item{package}{
    The package where the List class (\code{type}) is defined.
  }
  \item{regex}{
    A regular expression that uniquely identifies the file extension. 
  }
  \item{fnames}{
    A \code{character} vector of file names. 
  }
  \item{\dots}{
    Additional arguments passed to the List-class constructor
    (e.g., yieldSize for BamFileList).
  }
}

\details{

  \itemize{
    \item{registerFileType}{
      The \code{registerFileType} function adds entries to the file 
      type register created at load time. The point of the register
      is for discovery of file type (class) by file extension.
      These are List-type classes (e.g., BamFileList) that occupy
      the \code{fileList} slot of a GenomicFiles class.

      Each List class entry in the register is associated with
      (1) a regular expression that identifies the file extension,
      (2) a class and (3) the package where the class is defined.
      At load time the register is populated with classes known to
      GenomicFiles. New classes / file types can be added to
      the register with \code{registerFileType} by providing
      these three pieces of information.
    }
    \item{findTypeRegistry}{
      Searches the registry for a match to the extension of
      \code{fname}. Internal use only. 
    }
    \item{makeFileType}{
      Performs a look-up in the file registry based on the supplied
      regular expression; returns an object of the associated class. 
      Internal use only.
    }
  }
}

\value{
  \code{registerFileType}: NULL 
 
  \code{findTypeRegistry}: The regular expression associated with the file.

  \code{makeFileType}: A List-type object defined in the registry. 
}

\examples{

## At load time the registry is populated with file types
## known to GenomicFiles.
sapply(as.list(.fileTypeRegistry), "[", "type")

## Add a new class to the file register.
\dontrun{registerFileType(NewClassList, NewPackage, "\\.NewExtension$")}
}

\keyword{methods}
