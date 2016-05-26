\name{VcfStack}
\docType{class}

% Class:
\alias{class:VcfStack}
\alias{VcfStack-class}

\alias{Rangedclass:VcfStack}
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
\alias{seqinfo,VcfStack-method}
\alias{seqinfo<-,VcfStack-method}
\alias{colData,VcfStack-method}
\alias{colData<-,VcfStack,DataFrame-method}

\alias{rowRanges,RangedVcfStack-method}
\alias{rowRanges<-,RangedVcfStack,GRanges-method}

% Methods:
\alias{assay,VcfStack,missing-method}
\alias{assay,RangedVcfStack,missing-method}
\alias{readVcfStack}
\alias{show,VcfStack-method}

% Subsetting:

\alias{[,VcfStack,ANY,ANY-method}
\alias{[,VcfStack,ANY,ANY,ANY-method}

\alias{[,RangedVcfStack,missing,missing,missing-method}
\alias{[,RangedVcfStack,missing,character,missing-method}

% Helpers:
\alias{getVCFPath}
\alias{paths1kg}
\alias{dim,VcfStack-method}


\title{VcfStack and RangedVcfStack Objects}

\description{
  The \code{VcfStack} class is a vector of related VCF files, for instance 
  each file representing a seperate chromosome. The class helps manage these
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
	    A character vector of files paths pointing to VCF files.
	  }
	  \item{\code{seqinfo}}{
	    A \link[GenomeInfoDb]{Seqinfo} object describing the levels
            genome and circularity of each sequence.
	  }
	  \item{\code{colData}}{
	    An optional \link[S4Vectors]{DataFrame} describing each sample 
	    in the VcfStack.
	  }
	}
    }
    \item{}{
	\code{RangedVcfStack(vs=VcfStack(), rowRanges = GRanges())}
	Creates a RangedVcfStack object.
	\describe{
	  \item{\code{vs}}{
	    A \code{VcfStack} object.
	  }
	  \item{\code{rowRanges}}{
	    An optional \link[GenomicRanges]{GRanges} object associating the
	    genomic ranges of interest to the collection of VCF files.
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
      Get the short names of the files in \code{VcfStack}.
    }
    \item{dimnames(x))}{
      Get the sample names and the short names of files in \code{VcfStack}.
    }
    \item{files(x, \dots), files(x, \dots) <- value}{
      Get or set the files on \code{x}. \code{value} can be a character() of
      file paths or a \link[VariantAnnotation]{VcfFileList}.
    }
    \item{seqinfo(x), seqinfo(x, new2old = NULL, force = FALSE) <- value}{
      Get or set the seqinfo on \code{x}. See \link[GenomeInfoDb]{seqinfo<-} 
      for details on \code{new2old} and \code{force}. 
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
  In the code below, \code{x} is a VcfStack or RangedVcfStack object.

  \describe{
%    \item{assay(x, i, \dots)}{
%      Get matrix of genotype calls from the VCF files. 
%      See \link[VariantAnnotation]{genotypeToSnpMatrix}. Argument \code{i} 
%      must be missing. 
%    }
    \item{readVcfStack(x, i, j=colnames(x))}{
      Get content of VCF files in the VcfStack. \code{i} can be missing (if
      \code{x} is a \code{RangedVcfStack}) or a \code{GRanges} object. \code{j}
      can be missing or a character() vector of sample names (see
      \link[VariantAnnotation:VCFHeader-class]{samples}) present in the VCF
      files. 
    }
    \item{show(object)}{
      Display abbreviated information about \code{VcfStack} or 
      \code{RangedVcfStack} object.
    }
  }
}

\section{Subsetting}{
  In the code below, \code{x} is a VcfStack or RangedVcfStack object.

  \describe{
    \item{x[i, j]}{
      Get elements from ranges \code{i} and samples \code{j} as a
      \link[VariantAnnotation]{VCF} object. \code{i} can be missing or a 
      \code{GRanges} object. \code{j} can be missing or a character() vector of 
      sample names. 
    }    
  }
}

\section{Helpers}{
  
  \describe{
    \item{getVCFPath(vs, chrtok)}{
     
    }
    \item{paths1kg(chrtoks)}{
      Translate seqnames \code{chrtoks} to 1000 genomes genotypes vcf urls. 
    }
  }
}


\seealso{
  \link[VariantAnnotation]{VcfFile}, \link[VariantAnnotation]{VcfFileList}.
}

\author{Lori Shepherd {\url{mailto:Lori.Shepherd@RoswellPark.org}} and 
        Martin Morgan {\url{mailto:Martin.Morgan@RoswellPark.org}}}

\examples{
## Construction 
##
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

## Accessors
dim(stack)
colnames(stack)
rownames(stack)
dimnames(stack)
files(stack)
seqinfo(stack)
colData(stack)

## Methods
readVcfStack(stack, i=GRanges("20:862167-62858306"))
i = as(c("20:862167-62858306", "7:1-159138000"), "GRanges")
readVcfStack(stack, i=i, j="NA12891")

## Subsetting
stack[4:6,1:2]
stack[c("7","11","X"),"NA12891"]
stack[c("7","11","X"),2:3]

stack[GRanges("20:862167-62858306")]

## Helpers
getVCFPath(stack, 1)
paths1kg(1:3)
}


