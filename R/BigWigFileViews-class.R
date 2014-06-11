### =========================================================================
### BigWigFileViews methods
### =========================================================================

setClass("BigWigFileViews", contains="GenomicFileViews")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

setMethod(.validity, "BigWigFileViews", 
    function(object) {
        msg <- NULL
        if (!is(fileList(object), "BigWigFileList"))
            msg <- "'fileList' must be a 'BigWigFileList'"
        if (is.null(msg)) TRUE else msg
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generic 
###

setGeneric("BigWigFileViews",
           function(fileList,
                    fileSample=DataFrame(row.names=
                      make.unique(basename(fileList))),
                    fileRange=GRanges(),
                    fileExperiment=list(), 
                    yieldSize="NA_integer_",
                   .views_on_file="environment", ...)
           standardGeneric("BigWigFileViews"),
           signature="fileList")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

setMethod(BigWigFileViews, "missing", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(path(fileList)))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    stop("'fileList' must be character() or BigWigFileList")
})

setMethod(BigWigFileViews, "BigWigFileList", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(path(fileList)))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    new("BigWigFileViews", ..., 
        fileList=fileList,
        fileSample=fileSample, fileRange=fileRange,
        fileExperiment=fileExperiment, yieldSize=yieldSize, 
        .views_on_file=.views_on_file)
})

setMethod(BigWigFileViews, "character", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(fileList))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    new("BigWigFileViews", ..., 
        fileList=BigWigFileList(fileList),
        fileSample=fileSample, fileRange=fileRange,
        fileExperiment=fileExperiment, yieldSize=yieldSize, 
        .views_on_file=.views_on_file)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### base::mean is an S3 generic
###

## Needed for summary(..., type=mean)
setGeneric("mean", signature="x")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### coverage() and summary() methods
###

setMethod(coverage, "BigWigFileViews",
    function(x, ..., by="file", summarize=TRUE, as="RleList")
{
    MAPPER <- function(RANGE, FILE, ...) {
        if (as == "RleList") 
            import(FILE, which=RANGE, as=as)[RANGE]
        else 
            import(FILE, which=RANGE, as=as)
    }

    ## NOTE: Data in 'assays' of an SE object must be
    ##       organized as ranges by files. Currently
    ##       'assays' can hold a matrix of list elements
    ##       in the proper dimensions of ranges x files.
    ##       Specifying 'by=range' requires refactoring of
    ##       the list to fit in 'assays'. For this simple
    ##       coverage method 'by=range' and 'by=file' give the
    ##       same results. Choose 'by=file' which is faster.
    if (summarize) { 
        REDUCER <- function(mapped, ...) { mapped } 
        .summarizeView(x, MAPPER, REDUCER, ..., BY="file", as=as)
    } else {
        REDUCER <- function(mapped, ...) do.call(c, mapped)
        .reduce(fileRange(x), fileList(x), MAPPER, REDUCER, 
                ..., BY=by, as=as)
    }
})

setMethod(summary, "BigWigFileViews",
    function(object, ..., by="file", summarize=TRUE) 
{
    MAPPER <- function(range, file, ...) {
        sumres <- summary(file, which=range, asRangedData=TRUE, ...)
        do.call(c, sumres)$score
    }
    REDUCER <- function(mapped, ...) do.call(c, mapped)

    if (summarize) 
        .summarizeView(object, MAPPER, REDUCER, ..., BY=by)
    else 
        .reduce(fileRange(object), fileList(object), 
                MAPPER, REDUCER, ..., BY=by)
})
