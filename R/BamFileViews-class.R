### =========================================================================
### BamFileViews methods
### =========================================================================

setClass("BamFileViews", contains="GenomicFileViews")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generic 
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


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

.msg_BFV <- paste0("'BamFileViews' objects are defunct. ",
                   "Use 'GenomicFiles()' instead.")

setMethod(BamFileViews, "ANY", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(path(fileList)))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    .Defunct(msg=.msg_BFV)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### summarizeOverlaps(), countBam(), scanBam() 
###


setMethod(scanBam, "BamFileViews",
    function(file, index=file, ..., param=ScanBamParam(what=scanBamWhat()))
        .Defunct(msg=.msg_BFV)
)

setMethod(countBam, "BamFileViews",
    function(file, index=file, ..., param=ScanBamParam())
        .Defunct(msg=.msg_BFV)
)

setMethod("summarizeOverlaps", c("GRanges", "BamFileViews"),
    function(features, reads, mode=Union,
             ignore.strand=FALSE, ...)
        .Defunct(msg=.msg_BFV)
)

setMethod("summarizeOverlaps", c("GRangesList", "BamFileViews"),
    function(features, reads, mode=Union,
             ignore.strand=FALSE, ...)
        .Defunct(msg=.msg_BFV)
)
