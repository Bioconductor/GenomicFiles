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
### Validity 
###

setMethod(.validity, "GenomicFileViews", 
    function(object) 
{
        msg <- NULL
        if (length(fileList(object)) != nrow(fileSample(object)))
            msg <- c(msg,
            "length(fileList(object)) != nrow(fileSample(object))")
        if (is.null(msg)) TRUE else msg
})

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

setMethod("[", c("GenomicFileViews", "ANY", "missing"),
    function(x, i, j, ..., drop=TRUE)
        initialize(x, fileRange=fileRange(x)[i,])
)

setMethod("[", c("GenomicFileViews", "missing", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (is.character(j))
        j <- match(j, colnames(x))
    if (any(is.na(j)))
        stop("subscript 'j' out of bounds")
    initialize(x, 
               fileList=fileList(x)[j],
               fileSample=fileSample(x)[j,,drop=FALSE])
})

setMethod("[", c("GenomicFileViews", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (is.character(i))
        j <- match(i, rownames(x))
    if (is.character(j))
        j <- match(j, colnames(x))
    if (any(is.na(i)))
        stop("subscript 'i' out of bounds")
    if (any(is.na(j)))
        stop("subscript 'j' out of bounds")
    initialize(x, 
        fileRange=fileRange(x)[i,],
        fileList=fileList(x)[j],
        fileSample=fileSample(x)[j,,drop=FALSE])
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod(show, "GenomicFileViews", 
    function(object) 
{
    cat(class(object), "dim:",
        paste(dim(object), c("ranges", "samples"), collapse=" x "),
        "\n")
    cat("names:", BiocGenerics:::selectSome(names(object)), "\n")
    cat("detail: use fileList(), fileSample(), fileRange(), ...",
        "\n")
})
