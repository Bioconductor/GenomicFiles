### =========================================================================
### VcfStack and RangedVcfStack class
### =========================================================================

.validVcfStack = function(object) 
{
    msg <- NULL
    pa <- files(object)
    npa <- names(pa)
    if (!is.character(npa) || length(npa)!=length(pa))
        msg <- c(msg, "'files' must be character with one element per path")
    ## FIXME: enforce @files %in% seqnames(@rowRanges)?
    ## FIXME: enforce length(@files) == length(@seqinfo)?
    ## FIXME: colData requirements?
    if (is.null(msg)) TRUE else msg
}

setClass("VcfStack", 
    representation(
        files="character",
        seqinfo="Seqinfo",
        colData="DataFrame"
    ),
    validity=.validVcfStack  # verify 1-1 mapping from files to seqinfo
)

setClass("RangedVcfStack", 
    contains="VcfStack",
    representation(
        rowRanges="GRanges", 
        sampleNames="character"
    )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

VcfStack <- function(files, seqinfo, colData=DataFrame(), set.seqlstyle="NCBI") 
{
    if (!is(seqinfo, "Seqinfo")) 
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    pt <- files
    sn <- names(files)
    si <- seqinfo[sn]
    seqlevelsStyle(si) <- set.seqlstyle
    names(files) <- seqnames(si)
    tmp <- new("VcfStack", files=files, colData=colData, seqinfo=si)
    tmp
}

RangedVcfStack <- function(vs, rowRanges = GRanges(), sampleNames=character()) 
{
    stopifnot(is(vs, "VcfStack"))
    new("RangedVcfStack", vs, rowRanges=rowRanges, sampleNames=sampleNames)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters 
###

setMethod("files", "VcfStack",
    function(x, ...) x@files 
)
setReplaceMethod("files", c("VcfStack", "character"),
    function(x, ..., value)
{
    x@files <- value 
    x
})

## seqinfo (also seqlevels, genome, seqlevels<-, genome<-), seqinfo<-
setMethod(seqinfo, "VcfStack",
    function(x) x@seqinfo 
)
setReplaceMethod("seqinfo", "VcfStack",
    function (x, new2old = NULL, force = FALSE, value)
{
    if (!is(value, "Seqinfo")) 
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    x@seqinfo <- value
    x
})

setMethod(colData, "VcfStack",
    function(x) x@colData 
)
setReplaceMethod("colData", c("VcfStack", "DataFrame"),
    function(x, ..., value)
{
    x@colData <- value 
    x
})

## Replaces "bindRanges" and "bindRanges<-"
setMethod("rowRanges", "RangedVcfStack", 
    function(x, ...) x@rowRanges
)
setReplaceMethod("rowRanges", c("RangedVcfStack", "GRanges"), 
    function(x, ..., value) 
{
        x@rowRanges <- value 
        x 
})

## Replaces "bindSampleNames" and "bindSampleNames<-"
setMethod("sampleNames", "RangedVcfStack", 
    function(object) 
{
        sn <- object@sampleNames 
        if (length(sn)) 
            sn
        else
            samples(scanVcfHeader(files(object)[[1]]))
})
setReplaceMethod("sampleNames", c("RangedVcfStack", "character"), 
    function(object, value) 
{
        object@sampleNames <- value 
        object 
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods
###

setGeneric("subsetSample", function(x, j, ...) standardGeneric("subsetSample"))
setMethod("subsetSample", c("RangedVcfStack"), 
    function(x, j, ...) 
{
    sampleNames(x)  <- intersect(j, sampleNames(x))
    x
})

## FIXME: enforce 'i' numeric?
setMethod("assay", c("RangedVcfStack", "missing"), 
    function(x, i, ...) 
{
        sn <- sampleNames(x)
        rr <- rowRanges(x)
        if (length(sn))
            rr <- rr[sn] 
        mat <- readGT(files(x)[unique(seqnames(rr))], 
                      param=ScanVcfParam(which=rr), ...) 
        t(as(genotypeToSnpMatrix(mat)$genotypes, "numeric"))
})

readVcfStack <- function(x, i, j=character())
{
    stopifnot(is(x, "VcfStack"))
    stopifnot(is(i, "GenomicRanges"))
    qseqnames <- as.character(unique(seqnames(i)))
    if (length(qseqnames) != 1L)
        stop("seqnames(i) must have one unique seqname")

    ## FIXME: assume seqnames will match file paths?
    path2use <- files(x)[qseqnames]
    param <- ScanVcfParam(which=i, samples=j)
    readVcf(path2use, param=param, genome=genome(i)[1])
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###
 
setMethod("[", c("VcfStack", "ANY", "ANY", "missing"),
    function(x, i, j, drop)
{
        if (!missing(i)) {
            if (!is(i, "character"))
                stop("'i' must be a character vector of seqnames")
            if (!all(i %in% seqnames(x)))
                stop("seqnames in 'i' must be present in seqnames(x)")
            seqnames(x) <- seqinfo(x)[i] 
            ## FIXME: assume seqnames will match file paths?
            files(x) <- files(x)[i] 
        }
        if (!missing(j)) {
            ## FIXME: subset colData?
        }
        x
})

setMethod("[", c("RangedVcfStack", "ANY", "ANY", "missing"),
    function(x, i, j, drop) 
{
        if (!missing(i)) {
            if (!is(i, "character"))
                stop("'i' must be a character vector of seqnames")
            rowRanges(x) <- rowRanges(x)[i]
        }
        if (!missing(j)) {
            if (!is(j, "character"))
                stop("'j' must be a character vector of sample names")
            sampleNames(x) <- unique(j)
        }
        callNextMethod()
})

setMethod("[", c("RangedVcfStack", "missing", "missing", "missing"),
    function(x, i, j, drop) x
)  # above needed for biocMultiAssay validity method checker which runs x[]

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show() 
###

setMethod("show", "VcfStack", function(object) {
    cat("VcfStack instance with", length(files(object)), "files.\n")
    cat("Genome build recorded as ", genome(seqinfo(object))[1], ".\n", sep="")
    cat("use 'readVcfStack() to extract VariantAnnotation VCF.\n")
    show(seqinfo(object))
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers
###

getVCFPath <- function(vs, chrtok) {
    stopifnot(is.atomic(chrtok), length(chrtok)==1)
    files(vs)[chrtok]
}

paths1kg <- function(chrtoks) sapply(chrtoks, .path1kg, USE.NAMES=FALSE)
.path1kg <- function (chrtok) 
{
    stopifnot(length(chrtok)==1 && is.atomic(chrtok))
    if (is.numeric(chrtok)) 
        chrtok = as.integer(chrtok)
    if (is(chrtok, "integer")) 
        chrtok = paste0("chr", chrtok)
    if (length(grep("chr", chrtok)) < 1) 
        warning("probably need 'chr' in input string")
    tmplate = "http://1000genomes.s3.amazonaws.com/release/20130502/ALL.%%N%%.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
    if (length(grep("X", chrtok)) > 0) 
        tmplate = "http://1000genomes.s3.amazonaws.com/release/20130502/ALL.%%N%%.phase3_shapeit2_mvncall_integrated_v1b.20130502.genotypes.vcf.gz"
    if (length(grep("Y", chrtok)) > 0) 
        tmplate = "http://1000genomes.s3.amazonaws.com/release/20130502/ALL.%%N%%.phase3_integrated_v1b.20130502.genotypes.vcf.gz"
    ans = as.character(gsub("%%N%%", chrtok, tmplate))
    names(ans) = chrtok
    ans
}
