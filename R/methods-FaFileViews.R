### =========================================================================
### FaFileViews methods
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors
###

setMethod(FaFileViews, "missing", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(path(fileList)))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    stop("'fileList' must be character() or FaFileList")
})

setMethod(FaFileViews, "FaFileList", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(path(fileList)))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    new("FaFileViews", ...,
        fileList=fileList,
        fileSample=fileSample, fileRange=fileRange,
        fileExperiment=fileExperiment, yieldSize=yieldSize, 
        .views_on_file=.views_on_file)
})

setMethod(FaFileViews, "character", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(fileList))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    new("FaFileViews", ...,
        fileList=FaFileList(lapply(fileList, FaFile, ...)),
        fileSample=fileSample, fileRange=fileRange, 
        fileExperiment=fileExperiment, yieldSize=yieldSize, 
        .views_on_file=.views_on_file)
})
