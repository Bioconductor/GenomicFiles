### =========================================================================
### FaFileViews methods
### =========================================================================

setClass("FaFileViews", contains="GenomicFileViews")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generic 
###

setGeneric("FaFileViews",
           function(fileList,
                    fileSample=DataFrame(row.names=
                      make.unique(basename(fileList))),
                    fileRange=GRanges(),
                    fileExperiment=list(),
                    yieldSize="NA_integer_",
                   .views_on_file="environment", ...)
           standardGeneric("FaFileViews"),
           signature="fileList")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

setMethod(FaFileViews, "ANY", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(path(fileList)))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    .msg_FFV <- paste0("'FaFileViews()' is defunct. ",
                  "Use 'GenomicFiles()' instead.")
    .Defunct(msg=.msg_FFV)
})
