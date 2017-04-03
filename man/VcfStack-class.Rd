\name{VcfStack}
\docType{class}

% Class:
\alias{class:VcfStack}
\alias{VcfStack-class}

\alias{RangedVcfStack-class}

% Constructors:
\alias{VcfStack}
\alias{RangedVcfStack}

% Accessors:
\alias{colnames,VcfStack-method}
\alias{rownames,VcfStack-method}
\alias{dimnames,VcfStack-method}
\alias{files,VcfStack-method}
\alias{files<-,VcfStack,character-method}
\alias{files<-,VcfStack,VcfFile-method}
\alias{files<-,VcfStack,VcfFileList-method}
\alias{seqinfo,VcfStack-method}
\alias{seqinfo<-,VcfStack-method}
\alias{seqinfo<-,RangedVcfStack-method}
\alias{seqlevelsStyle<-,VcfStack-method}
\alias{seqlevelsStyle<-,RangedVcfStack-method}
\alias{colData,VcfStack-method}
\alias{colData<-,VcfStack,DataFrame-method}

\alias{rowRanges,RangedVcfStack-method}
\alias{rowRanges<-,RangedVcfStack,GRanges-method}

% Methods:
\alias{assay,VcfStack,ANY-method}
\alias{assay,RangedVcfStack,ANY-method}
\alias{readVcfStack}
\alias{show,VcfStack-method}

% Subsetting:
\alias{[,VcfStack,ANY,ANY-method}
\alias{[,VcfStack,ANY,ANY,ANY-method}
\alias{[,RangedVcfStack,ANY,ANY-method}
\alias{[,RangedVcfStack,ANY,ANY,ANY-method}

% Helpers:
\alias{paths1kg}
\alias{dim,VcfStack-method}

\title{VcfStack and RangedVcfStack Objects}

\description{
  The \code{VcfStack} class is a vector of related VCF files, for instance
  each file representing a separate chromosome. The class helps manage these
  files as a group. The \code{RangedVcfStack} class extends \code{VcfStack} by
  associating genomic ranges of interest to the collection of VCF files.
}

\section{Constructor}{
  \describe{
    \item{}{
        \code{VcfStack(files=NULL, seqinfo=NULL, colData=NULL)}
        Creates a VcfStack object.
        \describe{
          \item{\code{files}}{
	    A VcfFilelist object. If a VcfFile or character vector is
	    given a VcfFileList will be coerced. The character vector
	    should be files paths pointing to VCF files. The
            character vector must be named, with names correspond to
            seqnames in each VCF file.
          }
          \item{\code{seqinfo}}{
            A \link[GenomeInfoDb]{Seqinfo} object describing the levels
            genome and circularity of each sequence.
          }
          \item{\code{colData}}{
            An optional \link[S4Vectors]{DataFrame} describing each sample
            in the VcfStack. When present, row names must correspond to
            sample names in the VCF file.
          }
        }
    }
    \item{}{
        \code{RangedVcfStack(vs=NULL, rowRanges=NULL)}
        Creates a RangedVcfStack object.
        \describe{
          \item{\code{vs}}{
            A \code{VcfStack} object.
          }
          \item{\code{rowRanges}}{
            An optional \link[GenomicRanges]{GRanges} object associating
            the genomic ranges of interest to the collection of VCF
            files. The seqnames of \code{rowRanges} are a subset of
            \code{seqnames(vs)}. If missing, a default is created from
            the \code{seqinfo} object of the provided \code{VcfStack}.
          }
        }
    }
  }
}


\section{Accessors}{
  In the code below, \code{x} is a VcfStack or RangedVcfStack object.

  \describe{
    \item{dim(x)}{
      Get the number of files and samples in the \code{VcfStack} object.
    }
    \item{colnames(x, do.NULL=TRUE, prefix="col")}{
      Get the sample names in the \code{VcfStack}.
    }
    \item{rownames(x), do.NULL=TRUE, prefix="row")}{
      Get the names of the files in \code{VcfStack}.
    }
    \item{dimnames(x))}{
      Get the names of samples and the names of files in \code{VcfStack}.
    }
    \item{files(x, \dots), files(x, \dots) <- value}{
      Get or set the files on \code{x}. \code{value} can be a named
      character() of file paths or a
      \link[VariantAnnotation]{VcfFileList}. The return value will be a
      \link[VariantAnnotation]{VcfFileList}.
    }
    \item{seqinfo(x), seqinfo(x, new2old = NULL, force = FALSE) <- value}{
      Get or set the seqinfo on \code{x}. See \link[GenomeInfoDb]{seqinfo<-}
      for details on \code{new2old} and \code{force}.
    }
    \item{seqlevelsStyle(x) <- value}{
      Set the seqlevels according to the supplied style. File names and
      rowRanges will also be updated if applicable.
      See \link[GenomeInfoDb]{seqlevelsStyle<-} for more details.
    }    
    \item{colData(x), colData(x, \dots) <- value}{
      Get or set the \code{colData} on \code{x}. \code{value} is a
      \link[S4Vectors]{DataFrame}.
    }
    \item{rowRanges(x), rowRanges(x, \dots) <- value}{
      Get or set the \code{rowRanges} on \code{x}. \code{x} has to be a
      \code{RangedVcfStack} object. \code{value} is a
      \link[GenomicRanges]{GRanges}.
    }
  }
}

\section{Methods}{
  In the code below, \code{x} is a VcfStack or RangedVcfStack
  object. \code{i} is a \link[GenomicRanges]{GRanges} object,
  character() vector of \link[GenomeInfoDb:Seqinfo-class]{seqnames},
  numeric() vector, logical() vector, or can be missing. For a
  RangedVcfStack object, assay and readVcfStack will use the associated
  \code{rowRanges} object for \code{i}.
  \describe{
    \item{assay(x, i, \dots)}{
      Get matrix of genotype calls from the VCF files.
      See \link[VariantAnnotation]{genotypeToSnpMatrix}. Argument \code{i}
      specifies which files to read.
    }
    \item{readVcfStack(x, i, j=colnames(x), param=ScanVcfParam())}{
      Get content of VCF files in the VcfStack. \code{i} indicates which
      files to read. \code{j} can be missing or a character() vector of
      sample names (see
      \link[VariantAnnotation:VCFHeader-class]{samples}) present in the
      VCF files. \code{param} is a
      \link[VariantAnnotation]{ScanVcfParam} object. If \code{param} is
      used \code{i} and \code{j} are ignored. 
    }
    \item{show(object)}{
      Display abbreviated information about \code{VcfStack} or
      \code{RangedVcfStack} object.
    }
  }
}

\section{Subsetting}{
  In the code below, \code{x} is a VcfStack or RangedVcfStack
  object.

  \describe{
    \item{x[i, j]}{
      Get elements from ranges \code{i} and samples \code{j} as a
      VcfStack or RangedVcfStack object. Note: for a \code{RangedVcfStack},
      the \code{rowRanges} object will also be subset.

      \code{i} can be missing, a character() vector of
      \link[GenomeInfoDb:Seqinfo-class]{seqnames}, numeric() vector of
      indexes, logical() or \code{GRanges} object. When \code{i} is a
      \code{GRanges} object,  \code{seqnames(i)} is then
      used to subset the files in the VcfStack.

      \code{j} can be missing, a character() vector of sample names, a
      numeric(), logical() vector.

    }
  }
}

\section{Helpers}{

  \describe{
    \item{getVCFPath(vs, chrtok)}{
      Deprecated. Use \code{files(vs)[chrtok]} instead.
    }
    \item{paths1kg(chrtoks)}{
      Translate seqnames \code{chrtoks} to 1000 genomes genotype VCF urls.
    }
  }
}


\seealso{
  \link[VariantAnnotation]{VcfFile}, \link[VariantAnnotation]{VcfFileList}.
}

\author{Lori Shepherd {\url{mailto:Lori.Shepherd@RoswellPark.org}} and
        Martin Morgan {\url{mailto:Martin.Morgan@RoswellPark.org}}}

\examples{
## ---------------------------------------------------------------------
## CONSTRUCTION
## ---------------------------------------------------------------------
## point to VCF files and add names corresponding to the sequence
## present in the file
extdata <- system.file(package="GenomicFiles", "extdata")
files <- dir(extdata, pattern="^CEUtrio.*bgz$", full=TRUE)
names(files) <- sub(".*_([0-9XY]+).*", "\\\\1", basename(files))

## input data.frame describing the length of each sequence, coerce to
## 'Seqinfo' object
seqinfo <- as(readRDS(file.path(extdata, "seqinfo.rds")), "Seqinfo")

stack <- VcfStack(files, seqinfo)
stack

## Use seqinfo from VCF files instead of explict value
stack2 <- VcfStack(files)

rstack <- RangedVcfStack(stack)
gr <- GRanges(c("7:1-159138000", "X:1-155270560"))
rstack2 <- RangedVcfStack(stack, gr)
rstack2

## ---------------------------------------------------------------------
## ACCESSORS
## ---------------------------------------------------------------------
dim(stack)
colnames(stack)
rownames(stack)
dimnames(stack)
head(files(stack))
seqinfo(stack)
colData(stack)

## ---------------------------------------------------------------------
## METHODS
## ---------------------------------------------------------------------
readVcfStack(stack, i=GRanges("20:862167-62858306"))
i <- GRanges(c("20:862167-62858306", "7:1-159138000"))
readVcfStack(stack, i=i, j="NA12891")

head(assay(stack, gr))
head(assay(rstack2))

seqlevels(stack2)
rownames(stack2)
seqlevelsStyle(stack2)
seqlevelsStyle(stack2) <- "UCSC"
seqlevelsStyle(stack2)
seqlevels(stack2)
rownames(stack2)

## ---------------------------------------------------------------------
## SUBSETTING
## ---------------------------------------------------------------------
## select rows 4, 5, 6 and samples 1, 2
stack[4:6, 1:2]
## select rownames "7", "11" and sample "NA12891"
stack[c("7", "11"), "NA12891"]
stack[c("7", "11", "X"), 2:3]
## subset with GRanges
stack[GRanges("20:862167-62858306")]

rstack2[]
rstack2[,1]

## ---------------------------------------------------------------------
## HELPERS
## ---------------------------------------------------------------------
paths1kg(1:3)
}
