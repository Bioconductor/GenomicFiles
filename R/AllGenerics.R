### =========================================================================
### All generics 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

setGeneric("BamFileViews",
           function(fileList,
                    fileSample=DataFrame(row.names=
                      make.unique(basename(fileList))),
                    fileRange=GRanges(),
                    fileExperiment=list(), 
                    yieldSize="NA_integer_",
                   .views_on_file="environment", ...)
           standardGeneric("BamFileViews"),
           signature="fileList")

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

setGeneric("TabixFileViews",
           function(fileList,
                    fileSample=DataFrame(row.names=
                      make.unique(basename(fileList))),
                    fileRange=GRanges(),
                    fileExperiment=list(),
                    yieldSize="NA_integer_",
                   .views_on_file="environment", ...)
           standardGeneric("TabixFileViews"),
           signature="fileList")

setGeneric("VCFFileViews",
           function(fileList,
                    fileSample=DataFrame(row.names=
                      make.unique(basename(fileList))),
                    fileRange=GRanges(),
                    fileExperiment=list(),
                    yieldSize="NA_integer_",
                   .views_on_file="environment", ...)
           standardGeneric("VCFFileViews"),
           signature="fileList")

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
### reduceByFile and reduceByRange
###

setGeneric("reduceByRange", 
    function(ranges, files, MAPPER, REDUCER, iterate=FALSE, ..., init)
    standardGeneric("reduceByRange"),
    signature="ranges")

setGeneric("reduceByFile", 
    function(ranges, files, MAPPER, REDUCER, iterate=FALSE, ..., init)
    standardGeneric("reduceByFile"),
    signature="ranges")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### base::mean is an S3 generic
###

setGeneric("mean", signature="x")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### pack and unpack
###

setGeneric("pack", function(x, ...)
    standardGeneric("pack"),
    signature="x")

setGeneric("unpack", function(flesh, skeleton, ...)
    standardGeneric("unpack"),
    signature=c("flesh", "skeleton"))

