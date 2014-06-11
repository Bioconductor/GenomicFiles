### =========================================================================
### GenomicFiles class
### =========================================================================

setGeneric(".validity", function(object) standardGeneric(".validity"))

setClass("GenomicFiles",
    representation(
        files="ANY",
        rowData="DataFrame",
        colData="DataFrame",
        exptData="SimpleList"),
    prototype(
        files=character()),
    validity=.validity)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity 
###

setMethod(.validity, "GenomicFiles", 
    function(object) 
{
    msg <- NULL
    if (length(files(object)) != nrow(colData(object)))
        msg <- "'length(files(object))' must equal 'nrow(colData(object))'"

    msg 
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generic 
###

setGeneric("GenomicFiles",
    function(files, rowData=DataFrame(), colData=DataFrame(),
             exptData="SimpleList", ...) 
        standardGeneric("GenomicFiles"),
    signature="files")

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

setMethod(GenomicFiles, "character",
   function(files, rowData=DataFrame(), colData=DataFrame(),
            exptData=SimpleList(), ...)
{
    if (length(files)) {
        if (missing(colData))
            colData <- DataFrame(row.names=files)
        else
            rownames(colData) <- files
    }
    new("GenomicFiles", files=files, rowData=rowData, 
        colData=colData, exptData=exptData, ...)
})

setMethod(GenomicFiles, "List",
   function(files, rowData=DataFrame(), colData=DataFrame(),
            exptData=SimpleList(), ...)
{
    ## FIXME: require List be named?
    if (!is.null(nms <- names(files))) {
        if (missing(colData))
            colData <- DataFrame(row.names=nms)
        else
            rownames(colData) <- nms
    } 
    new("GenomicFiles", files=files, rowData=rowData, 
        colData=colData, exptData=exptData, ...)
})

setMethod(GenomicFiles, "list",
    function(files, ...)
{
    GenomicFiles(as(files, "List"), ...)
})

setMethod(GenomicFiles, "missing",
    function(files, ...)
{
    GenomicFiles(character(), ...)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Getters and Setters
###

setGeneric("files", function(x, ...) standardGeneric("files"))
setMethod("files", "GenomicFiles",
    function(x, ...)
{
    slot(x, "files")
})

setGeneric("files<-",
    function(x, ..., value) standardGeneric("files<-"))
setReplaceMethod("files", c("GenomicFiles", "ANY"),
    function(x, ..., value)
{
    if (!is(value, "character") && !is(value, "List") && !is(value, "list"))
        stop("'value' must be a 'character' or 'List'")
 
    if (is(value, "list"))
        value <- as(value, "List")
    if (is.null(nms <- names(value)))
        nms <- value
    rownames(colData(x)) <- nms 
    slot(x, "files") <- value
    x 
})

setMethod(rowData, "GenomicFiles",
    function(x, ...) slot(x, "rowData"))

setReplaceMethod("rowData", c("GenomicFiles", "DataFrame"),
    function(x, ..., value)
{
    slot(x, "rowData") <- value
    x
})

setMethod(colData, "GenomicFiles",
    function(x, ...) slot(x, "colData"))

setReplaceMethod("colData", c("GenomicFiles", "DataFrame"),
    function(x, ..., value)
{
    if (length(files(x)) != nrow(value))
        stop("'length(files(x))' must equal 'nrow(value)'")
    if (is.null(nms <- names(files(x))))
        nms <- files(x)
    rownames(value) <- nms 
    slot(x, "colData") <- value
    x
})

setMethod(exptData, "GenomicFiles",
    function(x, ...) slot(x, "exptData"))

setReplaceMethod("exptData", c("GenomicFiles", "SimpleList"),
    function(x, ..., value)
{
    slot(x, "exptData") <- value
    x 
})

setReplaceMethod("exptData", c("GenomicFiles", "list"),
    function(x, ..., value)
{
    slot(x, "exptData") <- SimpleList(value)
    x 
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### dim, dimnames, names
###

setMethod(dim, "GenomicFiles",
    function(x)
{
    c(nrow(rowData(x)), nrow(colData(x)))
})

setMethod(dimnames, "GenomicFiles",
    function(x)
{
    list(names(rowData(x)), rownames(colData(x)))
})

setReplaceMethod("dimnames", c("GenomicFiles", "list"),
    function(x, value)
{
    rowData <- rowData(x)
    names(rowData) <- value[[1]]
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    slot(x, "rowData") <- rowData
    slot(x, "colData") <- colData
})

setReplaceMethod("dimnames", c("GenomicFiles", "NULL"),
    function(x, value)
{
    dimnames(x) <- list(NULL, NULL)
    x
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

setMethod("[", c("GenomicFiles", "ANY", "missing"),
    function(x, i, j, ..., drop=TRUE)
        initialize(x, rowData=rowData(x)[i,])
)

setMethod("[", c("GenomicFiles", "missing", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (is.character(j))
        j <- match(j, colnames(x))
    if (any(is.na(j)))
        stop("subscript 'j' out of bounds")
    initialize(x, 
               rowData=rowData(x)[j],
               files=files(x)[j,,drop=FALSE])
})

setMethod("[", c("GenomicFiles", "ANY", "ANY"),
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

setMethod(show, "GenomicFiles", 
    function(object) 
{
    cat(class(object), "dim:",
        paste(dim(object), c("ranges", "files"), collapse=" x "),
        "\n")
    cat("names:", BiocGenerics:::selectSome(basename(files(object))), "\n")
    cat("detail: use files(), rowData(), colData(), ...",
        "\n")
})
