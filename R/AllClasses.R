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
setMethod(.validity, ".FileList", 
    function(object) {
        msg <- NULL
        if (!c("path", "index") %in% names(object))
            msg <- c(msg, paste0("names(fileList(object))",
                     " must include 'path' and 'index'"))
        if (is.null(msg)) TRUE else msg
    }
)
