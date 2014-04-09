### =========================================================================
### BigWigFileViews methods
### =========================================================================

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
### coverage() and summary() methods
###

setMethod(coverage, "BigWigFileViews",
    function(x, ..., by="file", summarize=TRUE, as="RleList")
{
    MAP <- function(FILE, RANGE, ...) {
        if (as == "RleList") 
            import(FILE, which=RANGE, as=as)[RANGE]
        else 
            import(FILE, which=RANGE, as=as)
    }
    REDUCE <- function(MAPPED, ...) do.call(c, MAPPED)

    ## NOTE: Data in 'assays' of an SE object must be
    ##       organized as ranges by files. Currently
    ##       'assays' can hold a 'list' where the
    ##       overall length matches the number of files
    ##       and elementLengths match the number of ranges.
    ##       Specifying 'by=range' requires refactoring of
    ##       the list to fit in 'assays'. For this simple
    ##       coverage method 'by=range' and 'by=file' give the
    ##       same results. Choose 'by=file' which is faster.
    if (summarize) 
        .summarizeView(x, MAP, REDUCE, ..., BY="file")
    else 
        .reduce(x, MAP, REDUCE, ..., BY=by)
})

setMethod(summary, "BigWigFileViews",
    function(object, ..., by="file", summarize=TRUE) 
{
    MAP <- function(FILE, RANGE, ...) {
        sumres <- summary(FILE, which=RANGE, asRangedData=TRUE, ...)
        do.call(c, sumres)$score
    }
    REDUCE <- function(MAPPED, ...) do.call(c, MAPPED)

    if (summarize) 
        .summarizeView(object, MAP, REDUCE, ..., BY=by)
    else 
        .reduce(object, MAP, REDUCE, ..., BY=by)
})
