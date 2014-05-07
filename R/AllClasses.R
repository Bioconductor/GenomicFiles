### =========================================================================
### All classes 
### =========================================================================

setGeneric(".validity", function(object) standardGeneric(".validity"))

setClass("GenomicFileViews",
    representation("VIRTUAL",
        fileList="List",
        fileSample="DataFrame",
        fileRange="GRanges",
        fileExperiment="list",
        yieldSize="integer",
        .views_on_file="environment"),
    prototype(
        fileList=List(),
        yieldSize=NA_integer_),
    validity=.validity)

setClass("BamFileViews", contains="GenomicFileViews")

setClass("FaFileViews", contains="GenomicFileViews")

setClass("TabixFileViews", contains="GenomicFileViews")

setClass("VcfFileViews", contains="GenomicFileViews")

setClass("BigWigFileViews", contains="GenomicFileViews")

.FileList <- setClass(".FileList", contains="SimpleList")

### FileManager
###
### Tasks to manage:
### - (split) divide problem for parallel eval
### - (apply / MAP) function executed on workers
### - (combine / REDUCE) collate results from workers 

### Thoughts:
### - Parallel exection:
###   Should be optional (large data, shared memory).
###
### - chunk/yield/iterate:
###   These are are independent of parallel and appropriate 
###   to include in apply/MAP step.
###
### - I/O considerations:
###   Potential problems with parallel reduceByRange() and 
###   single ff or hdf5 file (all hit same file, same time).
###
### - run byFile but reduce byRange:
###   Extraction of multiple ranges is more efficient
###   using a param (reduceByFile approach). If memory isn't 
###   a problem, we may consider a fast method to combine
###   elements of the list returned by bplapply() across
###   list elements. This would allow the user to run byFile
###   but reduce byRange.

setClass("FileManager",
    representation(
        files="character",
        split="ANY",        ## files, indices or ranges 
        apply="function",   ## pass extra args in ...
        combine="function", ## no extra args
        parallel="logical"),
    prototype(
        parallel=TRUE),
    validity=.validity)

### Sample constructor
### At minimum, user provides files and function. By default
### work is executed in parallel by file and unlisted.

##  FileManager(files,
##              apply,
##              split=files,
##              combine=CMBFUN <- function(x, ...) unlist(x),
##              parallel=TRUE, ...)
##
##  fm <- FileManager(c(file1, file2), countBam)
##  deploy(fm)
