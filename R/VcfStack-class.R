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
        colData="DataFrame",
        headers="list"
    ),
    validity=.validVcfStack  # verify 1-1 mapping from files to seqinfo
)

setClass("RangedVcfStack", 
    contains="VcfStack",
    representation(
        rowRanges="GRanges" 
    )
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

VcfStack <- function(files=character(), seqinfo=Seqinfo(), colData=DataFrame()) 
{
    if (!is(seqinfo, "Seqinfo")) 
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    pt <- files
    sn <- names(files)
    si <- seqinfo[sn]
    #seqlevelsStyle(si) <- set.seqlstyle
    names(files) <- seqnames(si)
    hs <- lapply(files, scanVcfHeader)
    tmp <- new("VcfStack", files=files, colData=colData, seqinfo=si,
       headers=hs)
    tmp
}

RangedVcfStack <- function(vs=VcfStack(), rowRanges = GRanges()) #, sampleNames=character()) 
{
    stopifnot(is(vs, "VcfStack"))
    new("RangedVcfStack", vs, rowRanges=rowRanges) #, sampleNames=sampleNames)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters 
###

setMethod("colnames", "VcfStack", function(x, do.NULL=TRUE, prefix="col") {
  rownames(colData(x))
})

setMethod("rownames", "VcfStack", function(x, do.NULL=TRUE, prefix="row") {
  names(files(x)) # ugly
})

setMethod("dim", "VcfStack", function(x) {
  c(length(x@files), length(vcfSamples(x@headers[[1]])))
})
 
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods
###


## FIXME: enforce 'i' numeric?
setMethod("assay", c("VcfStack", "missing"),
    function(x, i, ...)
{
        vcfob <- readVcf(files(x), genome=genome(x), ...)
        t(as(genotypeToSnpMatrix(vcfob)$genotypes, "numeric"))
})

## FIXME: enforce 'i' numeric?
setMethod("assay", c("RangedVcfStack", "missing"), 
    function(x, i, ...) 
{
        rr <- rowRanges(x)
        usn = unique(seqnames(rr))
#        mat <- readGT(files(x)[usn],  # used to work but now genotypeTo... needs ref
#                      param=ScanVcfParam(which=rr), ...) 
        vcfob <- readVcf(files(x)[unique(seqnames(rr))], genome=genome(x)[usn],
                      param=ScanVcfParam(which=rr), ...) 
        t(as(genotypeToSnpMatrix(vcfob)$genotypes, "numeric"))
})

readVcfStack <- function(x, i, j=character())
{
    stopifnot(is(x, "VcfStack"))
    if (is(x, "RangedVcfStack") && missing(i))
    i = rowRanges(x)
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
 
# we are going to allow 'row'-like subsetting using a GRanges
# sample subsetting is typical
# seqnames(i) picks the element of the VcfStack to read
# all subsetting is dependent on GRanges as row selection predicate
# if that is missing a message is given and the input returned

setMethod("[", c("VcfStack", "GenomicRanges", "character", "missing"),
   function(x, i, j, drop) {
    querseq = as.character(seqnames(i))
    stopifnot(length(unique(querseq))==1) # is this good enough?
    path2use = files(x)[querseq]
    param = ScanVcfParam(which=i)
    vcfSamples(param) = j
    readVcf(path2use, param=param, genome=genome(i)[1])
   })

setMethod("[", c("VcfStack", "GenomicRanges", "missing", "missing"),
   function(x, i, j, drop) {
    querseq = as.character(seqnames(i))
    stopifnot(length(unique(querseq))==1) # is this good enough?
    path2use = files(x)[querseq]
    param = ScanVcfParam(which=i)
    readVcf(path2use, param=param, genome=genome(i)[1])
   })

setMethod("[", c("VcfStack", "missing", "character", "missing"),
   function(x, i, j, drop) {
    message("omitting first subscript disallowed, please use GRanges subscripting")
    message("returning VcfStack unaltered.")
    x
})

setMethod("[", c("RangedVcfStack", "missing", "missing", "missing"),
    function(x, i, j, drop) {
#    message("omitting first subscript disallowed, please use GRanges subscripting")
#    message("returning VcfStack unaltered.")
#    x
#     i = rowRanges(x)
#     callNextMethod()
     x[ rowRanges(x), ]  # can't get callNext... dispatch right
})

setMethod("[", c("RangedVcfStack", "missing", "character", "missing"),
    function(x, i, j, drop) {
#    message("omitting first subscript disallowed, please use GRanges subscripting")
#    message("returning VcfStack unaltered.")
#    x
     i = rowRanges(x)
     callNextMethod()
})

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

