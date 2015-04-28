### =========================================================================
### GenomicFileViews (VIRTUAL)
### =========================================================================

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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters
###

fileList <- function(x, ...) 
    slot(x, "fileList")

`fileList<-` <-
    function(x, ..., value) initialize(x, fileList=value)

fileSample <-
    function(x, ...) slot(x, "fileSample")

`fileSample<-` <-
    function(x, ..., value) initialize(x, fileSample=value)
 
fileRange <-
    function(x, ...) slot(x, "fileRange")

`fileRange<-` <-
    function(x, ..., value) initialize(x, fileRange=value)

fileExperiment <-
    function(x, ...) slot(x, "fileExperiment")

`fileExperiment<-` <-
    function(x, ..., value) initialize(x, fileExperiment=value)

setMethod(yieldSize, "GenomicFileViews",
    function(object, ...) slot(object, "yieldSize"))

setReplaceMethod("yieldSize", "GenomicFileViews",
    function(object, ..., value)
{
    if (1L != length(value))
        stop("'value' must be length 1")
    slot(object, "yieldSize") <- as.integer(value)
    object
})

setMethod(names, "GenomicFileViews", 
    function(x) rownames(fileSample(x))
)

setReplaceMethod("names", "GenomicFileViews", 
    function(x, value) 
{ 
    rownames(fileSample(x)) <- value
    x
})

setMethod(dim, "GenomicFileViews", 
    function(x) c(length(fileRange(x)), length(fileList(x)))
)

setMethod(dimnames, "GenomicFileViews", 
    function(x) 
        list(names(fileRange(x)), rownames(fileSample(x)))
)

setReplaceMethod("dimnames", "GenomicFileViews", 
    function(x, value) 
{
        names(fileRange(x)) <- value[[1]]
        rownames(fileSample(x)) <- value[[2]]
        x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

.msg_GFV <- paste0("'GenomicFileViews' objects are defunct. ",
              "Use 'GenomicFiles()' instead.")

setMethod("[", c("GenomicFileViews", "ANY", "missing"),
    function(x, i, j, ..., drop=TRUE)
        .Defunct(msg=.msg_GFV)
)

setMethod("[", c("GenomicFileViews", "missing", "ANY"),
    function(x, i, j, ..., drop=TRUE)
        .Defunct(msg=.msg_GFV)
)

setMethod("[", c("GenomicFileViews", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
        .Defunct(msg=.msg_GFV)
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod(show, "GenomicFileViews", 
    function(object) 
        .Defunct(msg=.msg_GFV)
)
