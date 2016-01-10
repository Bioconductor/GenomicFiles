\name{VcfStack}
\alias{VcfStack}
\alias{VcfStack-class}
\alias{RangedVcfStack-class}
\alias{paths1kg}
\alias{bindRanges}
\alias{bindSampleNames}
\alias{[,RangedVcfStack,GRanges,character,missing-method}
\alias{[,RangedVcfStack,missing,missing,missing-method}
\alias{[,RangedVcfStack,GRanges,ANY,missing-method}
\alias{[,VcfStack,GenomicRanges,character,missing-method}
\alias{[,VcfStack,GenomicRanges,missing,missing-method}
\alias{assay,RangedVcfStack,missing-method}
\alias{bindRanges,RangedVcfStack,GRanges-method}
\alias{bindRanges,VcfStack,GRanges-method}
\alias{bindSampleNames,RangedVcfStack,character-method}
\alias{features}
\alias{features,RangedVcfStack-method}
\alias{samples,RangedVcfStack-method}
\alias{subsetSample}
\alias{subsetSample,RangedVcfStack-method}
\alias{show,VcfStack-method}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
manage a collection of VCF files that are split over chromosomes for a set of samples
}
\description{
manage a collection of VCF files that are split over chromosomes for a set of samples; includes utility for obtaining URLs for 1000 genomes VCF in S3
}
\usage{
VcfStack(paths, seqinfo, colData = DataFrame(), set.seqlstyle = "NCBI")
paths1kg(chrtoks)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{paths}{
Character vector of paths to tabix-indexed VCF; no checking of indexing.
Short names should be present.
}
  \item{seqinfo}{
\code{\link[GenomeInfoDb]{Seqinfo-class}} instance that comprehends
the sequences corresponding to \code{paths}
}
  \item{colData}{
optional information on samples in \code{\link[S4Vectors]{DataFrame}}
format, no validity checking of correspondence at this time, but
a manual check should be provided for
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
 \item{chrtoks}{vector of tokens to select chromosomes from
1000 genomes S3 bucket VCF}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
%\author{
%%  ~~who you are~~
%}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
pa = paths1kg(paste0("chr", c(21:22,"Y")))
library(Homo.sapiens)  # necessary?
ob = VcfStack(pa, seqinfo(Homo.sapiens))
myr = GRanges("22", IRanges(20e6,20.01e6))
\dontrun{
if (!is.null(nsl("www.r-project.org"))) { # test for connection
    ob[ myr, "NA12878" ]
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ models }
