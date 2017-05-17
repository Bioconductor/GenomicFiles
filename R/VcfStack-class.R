### =========================================================================
### VcfStack and RangedVcfStack class
### =========================================================================

.validVcfStack = function(object)
{
    msg <- NULL

    if (!all(rownames(object) %in% seqlevels(object)))
        msg <- c(msg, "all rownames(object) must be in seqlevels(object)")

    if (length(files(object))) {
        smps = samples(scanVcfHeader(files(object)[[1]]))
        if (!all(colnames(object) %in% smps))
            msg <- c(msg,
                     "all colnames(object) must be sample names in VCF 'files'")
        samplesOk <- sapply(files(object), function(file) {
            setequal(samples(scanVcfHeader(file)), smps)
        })
        if (!all(samplesOk))
            msg <- c(msg, "sample names are not consistent between VCF 'files'")
    }

    if (is.null(msg)) TRUE else msg
}

setClass("VcfStack",
    representation(
        files="VcfFileList",
        seqinfo="Seqinfo",
        colData="DataFrame"
    ),
    validity=.validVcfStack
)

.validRangedVcfStack = function(object)
{
    msg <- NULL

    if (!identical(seqinfo(rowRanges(object)), seqinfo(object)))
        msg <- c(msg,
                 "seqinfo() on rowRanges() differs from seqinfo() on object")

    if (!all(seqnames(rowRanges(object)) %in% rownames(object)))
        msg <- c(msg, "not all 'GRanges' seqnames are in VcfStack")

    if (is.null(msg)) TRUE else msg
}

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

VcfStack <- function(files=NULL, seqinfo=NULL, colData=NULL, index=TRUE)
{
    stopifnot(is.logical(index), length(index) == 1L, !is.na(index))

    if (is.null(files)) {
        files <- VcfFileList()
        header <- NULL
    } else {
        if (!is(files, "VcfFileList"))
            files = VcfFileList(files)
        if (index)
            files = indexVcf(files)
        header <- scanVcfHeader(files[[1]])
        
    }

    if (is.null(seqinfo)) {
        seqinfo <- if (length(files)) {
            seqinfo(files)
        } else Seqinfo()
    }

    if (is.null(colData) && length(files)) {
        colData <- DataFrame(row.names=samples(header))
    } else {
        colData <- as(colData, "DataFrame")
    }

    if (is.null(rownames(colData)) && length(files))
         stop("specify rownames in 'colData'")

    new("VcfStack", files=files, colData=colData, seqinfo=seqinfo)
}

RangedVcfStack <- function(vs=NULL, rowRanges=NULL)
{
    if (is.null(vs) && is.null(rowRanges)) {
        vs <- VcfStack()
        rowRanges <- GRanges()
    } else {
        stopifnot(is(vs, "VcfStack"))
        if (is.null(rowRanges)){
            rowRanges <- GRanges(seqinfo(vs))
            if (any(!seqnames(rowRanges) %in% rownames(vs)))
                rowRanges <- rowRanges[seqnames(rowRanges) %in% rownames(vs)]
        }
        new2old <- match(seqlevels(vs), seqlevels(rowRanges))
        seqinfo(rowRanges, new2old=new2old) <- seqinfo(vs)
    }

    new("RangedVcfStack", vs, rowRanges=rowRanges)
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
    files(x) <- VcfFileList(value)
    x
})

setReplaceMethod("files", c("VcfStack", "VcfFile"),
    function(x, ..., value)
{
    files(x) <- VcfFileList(value)
    x
})

setReplaceMethod("files", c("VcfStack", "VcfFileList"),
    function(x, ..., value)
{
    value <- indexVcf(value)
    initialize(x, files=value)
})

## seqinfo (also seqlevels, genome, seqlevels<-, genome<-)
setMethod(seqinfo, "VcfStack",
    function(x) x@seqinfo
)

setReplaceMethod("seqinfo", "VcfStack",
    function (x, new2old = NULL, pruning.mode = c("error", "coarse", "fine", "tidy"), value)
{
    initialize(x, seqinfo=value)
})

## H.P. 2017-04-29: I renamed 'force' -> 'pruning.mode'. Surprisingly this
## argument is ignored. That doesn't seem right.
setReplaceMethod("seqinfo", "RangedVcfStack",
    function (x, new2old = NULL, pruning.mode = c("error", "coarse", "fine", "tidy"), value)
{
    if (!is(value, "Seqinfo"))
        stop("the supplied 'seqinfo' must be a Seqinfo object")
    if (is.null(new2old))
        new2old <- match(seqnames(value), seqlevels(rowRanges(x)))
    rowRanges <- rowRanges(x)
    seqinfo(rowRanges, new2old=new2old) <- value
    initialize(x, seqinfo=value, rowRanges=rowRanges)
})

setReplaceMethod("seqlevelsStyle", "VcfStack",
    function(x, value)
{
    newSeqInfo <- seqinfo(x)
    seqlevelsStyle(newSeqInfo) <- value
    newFiles <- files(x)
    nms = names(newFiles)
    seqlevelsStyle(nms) <- value
    names(newFiles) <- nms
    initialize(x, seqinfo=newSeqInfo, files=newFiles)
})

setReplaceMethod("seqlevelsStyle", "RangedVcfStack",
    function(x, value)
{
    newSeqInfo <- seqinfo(x)
    seqlevelsStyle(newSeqInfo) <- value
    newFiles <- files(x)
    nms = names(newFiles)
    seqlevelsStyle(nms) <- value
    names(newFiles) <- nms
    newRange <- rowRanges(x)
    seqlevelsStyle(newRange) <- value
    initialize(x, seqinfo=newSeqInfo, files=newFiles, rowRanges=newRange)
})

setMethod(colData, "VcfStack",
    function(x) x@colData
)

setReplaceMethod("colData", c("VcfStack", "DataFrame"),
    function(x, ..., value)
{
    initialize(x, colData=value)
})

setMethod("rowRanges", "RangedVcfStack",
    function(x, ...) x@rowRanges
)

setReplaceMethod("rowRanges", c("RangedVcfStack", "GRanges"),
    function(x, ..., value)
{
    new2old <- match(seqlevels(x), seqlevels(value))
    seqinfo(value, new2old=new2old) <- seqinfo(x)
    initialize(x, rowRanges=value)
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

setMethod("assay", c("RangedVcfStack", "ANY"),
    function(x, i, ...)
{
    if (!missing(i))
        message(paste(strwrap(
            "RangedVcfStack uses rowRanges to subset; ignoring 'i'",
            exdent=4), collapse="\n"))
    i <- rowRanges(x)
    callNextMethod(x=x, i=i)
})

readVcfStack <- function(x, i, j=colnames(x), param=ScanVcfParam())
{
    stopifnot(is(x, "VcfStack"))
    if ((!missing(i) || !missing(j)) && !missing(param))
        stop("'i' and 'j' cannot be used with 'param'")

    gr <- NULL
    if (missing(param) && missing(i) && is(x, "RangedVcfStack")) {
        gr <- rowRanges(x)
    } else if (missing(param) && missing(i)) {
        gr <- GRanges()
    } else if (missing(param) && is(i, "GRanges")) {
        x <- x[unique(seqnames(i))]
        gr <- i
    } else if (missing(param)) {
        if (is.numeric(i))
            i = names(files(x))[i]
        x <- x[i]
        gr <- GRanges()
    } else {                            # use param
        gr <- GRanges(vcfWhich(param))
    }
    
    if (is.numeric(j)) {
        j <- colnames(x)[j]
    } else if (!missing(param)) {
        j <- vcfSamples(param)
    }

    genome <- genome(x)
    vcfSamples(param) <- j
    vcfWhich(param) <- gr

    idx <- if (length(gr) > 0){
        intersect(names(files(x)),
                  as.character(seqnames(gr))) } else {
                      names(files(x)) }
    
    vcf <- lapply(idx, function(i, files, genome, param) {
        file <- files[[i]]
        if (length(vcfWhich(param)))
            vcfWhich(param) <- vcfWhich(param)[i]
        readVcf(file, genome, param)
    }, files(x), genome, param)

    do.call(rbind, vcf)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

setMethod("[", c("VcfStack", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE){

    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,VcfStack,ANY,ANY-method'")

    if (missing(i) && missing(j)) {
        x
    } else if (missing(j)) {
        if (is(i, "GRanges")) {
            i <- as.character(seqnames(i))
        }
        initialize(x, files=files(x)[i])
    } else if (missing(i)) {
        initialize(x, colData=colData(x)[j,,drop=FALSE])
    } else {
        if (is(i, "GRanges")) {
            i <- as.character(seqnames(i))
        }
        initialize(x, files=files(x)[i], colData=colData(x)[j,,drop=FALSE])
    }
})

setMethod("[", c("RangedVcfStack", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE) {

    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,RangedVcfStack,ANY,ANY-method'")

    if (missing(i)) {
        i <- rownames(x)
    } else if (is(i, "GenomicRanges")) {
        stopifnot(all(seqnames(i) %in% rownames(x)))
        rowRanges(x) <- intersect(rowRanges(x), i)
    } else if (is(i, "character")) {
        stopifnot(all(i %in% rownames(x)))
        value <- rowRanges(x)
        rowRanges(x) <- value[seqnames(value) %in% i]
    } else {
        stopifnot(is(i, "numeric") || is(i, "logical"))
        value <- rowRanges(x)
        rowRanges(x) <- value[seqnames(value) %in% rownames(x)[i]]
    }

    if (missing(j))
        j <- colnames(x)

    callNextMethod(x=x, i=i, j=j)
})


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### show()
###

setMethod("show", "VcfStack", function(object) {
    cat("VcfStack object with ", nrow(object), " files and ",
        ncol(object), " samples",
        "\n", sep="")
    if (is(object, "RangedVcfStack")) {
        cat(summary(rowRanges(object)), "\n")
    }
    cat("Seqinfo object with", summary(seqinfo(object)), "\n")
    cat("use 'readVcfStack()' to extract VariantAnnotation VCF.\n")
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers
###

getVCFPath <- function(vs, chrtok) {
    .Deprecated("files(vs)[chrtok]")
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
