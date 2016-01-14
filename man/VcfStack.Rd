\name{VcfStack}

%% Classes
\alias{VcfStack}
\alias{VcfStack-class}
\alias{RangedVcfStack}
\alias{RangedVcfStack-class}

%% Subsetting
\alias{[,VcfStack,ANY,ANY,missing-method}
\alias{[,RangedVcfStack,ANY,ANY,missing-method}
\alias{[,RangedVcfStack,missing,missing,missing-method}

%% Accessors
\alias{files}
\alias{files,VcfStack-method}
\alias{files<-,VcfStack,character-method}
\alias{seqinfo}
\alias{seqinfo,VcfStack-method}
\alias{seqinfo<-,VcfStack-method}
\alias{colData}
\alias{colData,VcfStack-method}
\alias{colData<-,VcfStack,DataFrame-method}
\alias{sampleNames}
\alias{sampleNames,RangedVcfStack-method}
\alias{sampleNames<-,RangedVcfStack,character-method}
\alias{rowRanges}
\alias{rowRanges,RangedVcfStack-method}
\alias{rowRanges<-,RangedVcfStack,GRanges-method}

%% Methods
\alias{subsetSample}
\alias{subsetSample,RangedVcfStack-method}
\alias{assay}
\alias{assay,RangedVcfStack,missing-method}
\alias{show,VcfStack-method}

%% Other
\alias{paths1kg}
\alias{getVCFPath}
\alias{readVcfStack}

\title{
  Manage a collection of VCF files split over chromosomes for a set of samples
}

\description{
  Manage a collection of VCF files split over chromosomes for a set 
  of samples; includes utility for obtaining URLs for 1000 genomes VCF in S3.
}

\usage{
  ## Constructors
  VcfStack(files, seqinfo, colData = DataFrame(), set.seqlstyle = "NCBI")
  RangedVcfStack(vs, rowRanges = GRanges(), sampleNames = character())

  ## Accessors for VcfStack and RangedVcfStack
  files(x, ...)
  seqinfo(x)
  colData(x, ...)

  ## Accessors for RangedVcfStack only
  rowRanges(x, ...)
  sampleNames(object)

  ## Other
  assay(x, i, ...)
  subsetSample(x, j, ...)
  readVcfStack(x, i, j = character())
  getVCFPath(vs, chrtok)
  paths1kg(chrtoks)
}

\arguments{
  \item{files}{
    Character vector of paths to tabix-indexed VCF files; no checking of 
    indexing. Short names should be present.
  }
  \item{seqinfo}{
    \code{\link[GenomeInfoDb]{Seqinfo-class}} instance that comprehends
    the sequences corresponding to \code{paths}
  }
  \item{colData}{
    Optional information on samples in \code{\link[S4Vectors]{DataFrame}}
    format, no validity checking of correspondence at this time.
  }
  \item{vs}{
    \code{VcfStack} object.
  }
  \item{sampleNames}{
    Character vector of samples names to extract from VCF files.
  }
  \item{rowRanges}{
    \code{\link[GenomicRanges]{GRanges-class}} of regions to extract from 
    VCF files.
  }
  \item{set.seqlstyle}{
    Character atom, typically either \code{"NCBI"} or \code{"UCSC"},
    telling how to define \code{\link[GenomeInfoDb]{seqlevelsStyle}}
    for the \code{seqinfo} slot contents.  In the case of 1000 genomes,
    filenames have chr[nn], but sequence data within VCF lacks "chr"
    prefix.  Thus our queries must be in NCBI style, but our filename
    references include UCSC style.  See the example for one feasible
    approach.
  }
  \item{x}{
    \code{VcfStack} or \code{RangedVcfStack} object.
  }
  \item{object}{
    \code{VcfStack} or \code{RangedVcfStack} object.
  }
  \item{i}{
    Row-like index for \code{VcfStack}. For \code{readVcfStack} 
    \code{i} is a \code{GRanges}. Not currently used in \code{assay}.
  }
  \item{j}{
    Column-like or sample index for \code{VcfStack}. Usually a vector
    of sample names.
  }
  \item{chrtoks}{
    Vector of tokens to select chromosomes from 1000 genomes S3 bucket VCF.
  }
  \item{chrtok}{
    Single element of \code{chrtoks} vector.
  }
  \item{\dots}{
    Arguments passed to other methods.
  }
}

\examples{
pa <- paths1kg(paste0("chr", c(21:22,"Y")))
library(Homo.sapiens)  # necessary?
ob <- VcfStack(pa, seqinfo(Homo.sapiens))
  
\dontrun{ ## network connectivity
  if (!is.null(nsl("www.r-project.org"))) { # test for connection
    ## Read chromosome 22 file, all samples:
    ranges22 <- GRanges("22", IRanges(20e6,20.01e6))
    vcf22 <- readVcfStack(ob, ranges22)
    vcf22
    
    ## Read chromosome 21 file, sample NA12878:
    ranges21 <- GRanges("22", IRanges(20e6,20.01e6))
    vcf21_NA12878 <- readVcfStack(ob, ranges21, "NA12878")
    vcf21_NA12878
  }
}
}

\keyword{ models }
