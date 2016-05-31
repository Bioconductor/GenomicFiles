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

    if (any(!rownames(object) %in% seqlevels(object)))
        msg <- c(msg, "rownames must be in seqinfo object")

    if (length(files(object))) {
        smps = samples(scanVcfHeader(files(object)[1]))
        if (any(!colnames(object) %in% smps))
            msg <- c(msg, "a sample in 'colData' is not a sample in 'files'")
    
        if (!all(sapply(lapply(files(object),scanVcfHeader),samples) == smps))
            msg <- c(msg, "sample names are not consistent between 'files'")
    }

    if (is.null(msg)) TRUE else msg
}


.validRangedVcfStack = function(object)
{
    msg <- NULL
    if (any(!seqnames(rowRanges(object)) %in% rownames(object)))
        msg <- c(msg, "A 'GRange' seqname is not in VcfStack")

    if (is.null(msg)) TRUE else msg
}

setClass("VcfStack",
    representation(
        files="character",
        seqinfo="Seqinfo",
        colData="DataFrame"
    ),
    validity=.validVcfStack
)

setClass("RangedVcfStack",
    contains="VcfStack",
    representation(
        rowRanges="GRanges"
    ),
    validity=.validRangedVcfStack
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
##

VcfStack <- function(files=NULL, seqinfo=NULL, colData=NULL)
{
    if (is.null(files)) {
        files <- structure(character(), .Names=character())
        header <- NULL
    } else {
        header <- scanVcfHeader(files[1])
    }

    if (is.null(seqinfo)) {
        seqinfo <- if (length(files)) {
            seqinfo(header)
        } else Seqinfo()
    }

    pt <- files
    sn <- names(files)
    si <- seqinfo

    if (is.null(colData) && length(files)) {
        colData <- DataFrame(row.names=samples(header))
    } else {
        colData <- as(colData, "DataFrame")
    }

    if (is.null(rownames(colData)) && length(files))
         stop("Specify rownames in 'colData'")

    new("VcfStack", files=files, colData=colData, seqinfo=si)
}

RangedVcfStack <- function(vs=NULL, rowRanges = NULL) 
{
    if (is.null(vs) && is.null(rowRanges)) {
        new("RangedVcfStack", VcfStack(), rowRanges=GRanges()) 
    } else {
        stopifnot(is(vs, "VcfStack"))
        if (is.null(rowRanges)){
            rowRanges <- as(seqinfo(vs), "GRanges")
            if (any(!seqnames(rowRanges) %in% rownames(vs)))
                rowRanges <- rowRanges[seqnames(rowRanges) %in% rownames(vs)]
        }         
        new("RangedVcfStack", vs, rowRanges=rowRanges)
    }
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and setters
###

setMethod("dimnames", "VcfStack", function(x){
    list(names(files(x)), rownames(colData(x)))
})

setMethod("dim", "VcfStack", function(x) {
  c(length(files(x)), nrow(colData(x)))
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

## seqinfo (also seqlevels, genome, seqlevels<-, genome<-)
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

setMethod("rowRanges", "RangedVcfStack",
    function(x, ...) x@rowRanges
)

setReplaceMethod("rowRanges", c("RangedVcfStack", "GRanges"),
    function(x, ..., value)
{
        x@rowRanges <- value
        validObject(x)
        x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Other methods
###

setMethod("assay", c("VcfStack", "ANY"),
     function(x, i, ...)
{
    if (is(i, "GRanges")) {
        files <- files(x)[as.character(seqnames(i))]
    } else {
        files <- if (missing(i)) files(x) else files(x)[i]
        i <- GRanges(seqinfo(x))[names(files)]
    }

    i <- splitAsList(i, seq_along(i))
    genotypes <- Map(function(file, grange, genome) {
        ## FIXME: readGeno or other more efficient input?
        vcf <- readVcf(file, genome, grange)
        t(as(genotypeToSnpMatrix(vcf)$genotypes, "numeric"))
    }, files, i, MoreArgs=list(genome=genome(x)))

    do.call(rbind, genotypes)
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

readVcfStack <- function(x, i, j=colnames(x))
{
    stopifnot(is(x, "VcfStack"))
    if (is(x, "RangedVcfStack") && missing(i)){
        ranges = rowRanges(x)
        i = as.character(unique(seqnames(ranges)))
    }
    if (missing(i)) {
        i = names(files(x))
    }
    if (is(i, "GRanges")) {
        i = as.character(unique(seqnames(i)))
    } 

    if (is.numeric(j)) {
        j = colnames(x)[j]
    }

    path2use <- files(x)[i]

    files = lapply(path2use, readVcf, genome=genome(x), 
                   param=ScanVcfParam(samples=j))

    do.call(rbind,files)
   
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

setMethod("[", c("VcfStack", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE){

    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,VCF,ANY,ANY-method'")

    if (missing(i) && missing(j)) {
        x
    } else if (missing(j)) {
        if (is(i, "GRanges")) {
            i = as.character(seqnames(i))
        }
        initialize(x, files=files(x)[i])
    } else if (missing(i)) {
	initialize(x, colData=colData(x)[j,])
    } else {
        if (is(i, "GRanges")) {
            i = as.character(seqnames(i))
        }        
  	initialize(x, files=files(x)[i], colData=colData(x)[j,])
    }      
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
    cat("VcfStack object with ", nrow(object), " files and ", 
        ncol(object), " samples",
        "\ngenome: ", unique(genome(object)),
        "\n", sep="")
    cat("Seqinfo object with", summary(seqinfo(object)), "\n")
    
    ## FIXME: replace when GenomicRanges summary() implemented 
    if( is(object, "RangedVcfStack") ){
        x <- rowRanges(object)
        x_class <- class(x)
        x_len <- length(x)
        x_mcols <- mcols(x)
        x_nmc <- if (is.null(x_mcols)) 0L else ncol(x_mcols)
        cat(x_class, " object with ",
            x_len, " ", ifelse(x_len == 1L, "range", "ranges"),
            " and ",
            x_nmc, " metadata ", ifelse(x_nmc == 1L, "column", "columns"),
            "\n", sep="")
    }
    cat("use 'readVcfStack() to extract VariantAnnotation VCF.\n")

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
