### =========================================================================
### VcfStack and RangedVcfStack class
### =========================================================================

.paths = function(ob) ob@paths
.validVcfStack = function(object) 
{
    pa = .paths(object)
    npa = names(pa)
    if (!is.character(npa) || length(npa)!=length(pa)) 
        return("object@paths must be character with one element per path")
    TRUE
}

setClass("VcfStack", 
    representation(
        paths="character",
        seqinfo="Seqinfo",
        colData="DataFrame"
    ),
    validity=.validVcfStack  # verify 1-1 mapping from paths to seqinfo
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

VcfStack = function(paths, seqinfo, colData=DataFrame(), set.seqlstyle="NCBI") 
{
    pt <- paths
    sn <- names(paths)
    si <- seqinfo[sn]
    seqlevelsStyle(si) <- set.seqlstyle
    names(paths) <- seqnames(si)
    tmp <- new("VcfStack", paths=paths, colData=colData, seqinfo=si)
    tmp
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods
###

setGeneric("bindRanges", function(vs, gr, ...) standardGeneric("bindRanges"))
setMethod("bindRanges", c("VcfStack", "GRanges"), 
    function(vs, gr, ...) new("RangedVcfStack", vs, rowRanges=gr)
)
setMethod("bindRanges", c("RangedVcfStack", "GRanges"), 
    function(vs, gr, ...) {
        vs@rowRanges = gr
        vs
    }
)

setGeneric("bindSampleNames", function(vs, sn, ...) 
    standardGeneric("bindSampleNames")
)
setMethod("bindSampleNames", c("RangedVcfStack", "character"), 
    function(vs, sn, ...) {
        vs@sampleNames = sn
        vs
    }
)

setMethod("samples", "RangedVcfStack", 
    function(object) {
        if (length(sn <- object@sampleNames)>0) return(sn)
        samples(scanVcfHeader(object@paths[[1]]))
    }
)

setGeneric("subsetSample", function(x, j, ...) standardGeneric("subsetSample"))
setMethod("subsetSample", c("RangedVcfStack"), 
    function(x, j, ...) {
        if (length(sn <- x@sampleNames)>0) {
            x@sampleNames = intersect(j, x@sampleNames)
            return(x)
        }
        sn = samples(scanVcfHeader(x@paths[[1]]))
        x@sampleNames=intersect(j, sn)
        x
    }
)

setGeneric("features", function(x)standardGeneric("features"))
setMethod("features", "RangedVcfStack", 
    function(x) x@rowRanges
)

setMethod("assay", c("RangedVcfStack", "missing"), 
    function(x, i) {
        sn = x@sampleNames
        vp = ScanVcfParam(which=x@rowRanges)
        if (length(sn)>0) rd = x[ x@rowRanges, sn ]
        else rd = x[ x@rowRanges, ]
        mat = as( genotypeToSnpMatrix(rd)$genotypes, "numeric" )
        t(mat)
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

setMethod("[", c("VcfStack", "GenomicRanges", "character", "missing"),
    function(x, i, j, drop) {
        querseq = as.character(seqnames(i))
        stopifnot(length(unique(querseq))==1) # is this good enough?
        path2use = .paths(x)[querseq]
        param = ScanVcfParam(which=i)
        vcfSamples(param) = j
        readVcf(path2use, param=param, genome=genome(i)[1])
    }
)

setMethod("[", c("VcfStack", "GenomicRanges", "missing", "missing"),
    function(x, i, j, drop) {
        querseq = as.character(seqnames(i))
        stopifnot(length(unique(querseq))==1) # is this good enough?
        path2use = .paths(x)[querseq]
        param = ScanVcfParam(which=i)
        readVcf(path2use, param=param, genome=genome(i)[1])
    }
)

setMethod("[", c("RangedVcfStack", "GRanges", "ANY", "missing"),
    function(x, i, j, drop) {
        x@rowRanges=i
        if(!missing(j)) x@sampleNames=j
        x
    }
)  # just updates slots endomorphically

setMethod("[", c("RangedVcfStack", "missing", "missing", "missing"),
    function(x, i, j, drop) x
)  # above needed for biocMultiAssay validity method checker which runs x[]

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show() 
###

setMethod("show", "VcfStack", function(object) {
    cat("VcfStack instance with", length(object@paths), "paths.\n")
    cat("Genome build recorded as ", 
        as.character(genome(object@seqinfo)[1]), ".\n", sep="")
    cat("use '[ [GRanges], [sampleids] ]' to extract VariantAnnotation VCF.\n")
    show(object@seqinfo)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers
###

getVCFPath <- function(vs, chrtok) {
    stopifnot(is.atomic(chrtok), length(chrtok)==1)
    .paths(vs)[chrtok]
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
