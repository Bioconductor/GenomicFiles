### =========================================================================
### BigWigFileViews methods
### =========================================================================

setClass("BigWigFileViews", contains="GenomicFileViews")

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

.msg_BWFV <- paste0("'BigWigFileViews' objects are defunct. ",
                    "Use 'GenomicFiles()' instead.")

setMethod(BigWigFileViews, "missing", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(path(fileList)))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    .Defunct(msg=.msg_BWFV)
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
    .Defunct(msg=.msg_BWFV)
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
    .Defunct(msg=.msg_BWFV)
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
        .Defunct(msg=.msg_BWFV)
)

setMethod(summary, "BigWigFileViews",
    function(object, ..., by="file", summarize=TRUE) 
        .Defunct(msg=.msg_BWFV)
)
