### =========================================================================
### GenomicFiles class
### =========================================================================

setGeneric(".validity", function(object) standardGeneric(".validity"))

setClass("GenomicFiles",
    representation(
        files="ANY",
        rowData="GenomicRangesORGRangesList",
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
    function(rowData, files, ...) 
        standardGeneric("GenomicFiles"),
    signature=c("rowData", "files"))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

setMethod(GenomicFiles, c("GenomicRangesORGRangesList", "character"),
   function(rowData, files, colData=DataFrame(), exptData=SimpleList(), ...)
{
    if (length(files)) {
        if (is.null(nms <- names(files)))
            nms <- basename(files)
        if (missing(colData))
            colData <- DataFrame(row.names=nms)
        else
            rownames(colData) <- nms 
    }
    new("GenomicFiles", rowData=rowData, files=files, 
        colData=colData, exptData=exptData)
})

setMethod(GenomicFiles, c("GenomicRangesORGRangesList", "List"),
   function(rowData, files, colData=DataFrame(), exptData=SimpleList(), ...)
{
    if (length(files)) {
        if (is.null(nms <- names(files)))
            stop("'List' of files must be named")

        if (missing(colData))
            colData <- DataFrame(row.names=basename(nms))
        else
            rownames(colData) <- basename(nms)
    }
    new("GenomicFiles", rowData=rowData, files=files, 
        colData=colData, exptData=exptData)
})

setMethod(GenomicFiles, c("GenomicRangesORGRangesList", "list"),
    function(rowData, files, ...)
{
    GenomicFiles(rowData, as(files, "List"), ...)
})

setMethod(GenomicFiles, c("missing", "ANY"),
    function(rowData, files, ...)
{
    GenomicFiles(GRanges(), files, ...)
})

setMethod(GenomicFiles, c("missing", "missing"),
    function(rowData, files, ...)
{
    GenomicFiles(DataFrame(), character(), ...)
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
setReplaceMethod("files", c("GenomicFiles", "character"),
    function(x, ..., value)
{
    if (is.null(nms <- names(value)))
        nms <- value
    rownames(colData(x)) <- nms 
    slot(x, "files") <- value
    x 
})

setReplaceMethod("files", c("GenomicFiles", "List"),
    function(x, ..., value)
{
    if (is.null(nms <- names(value)))
        nms <- value
    rownames(colData(x)) <- nms 
    slot(x, "files") <- value
    x 
})

setMethod(rowRanges, "GenomicFiles",
    function(x, ...) slot(x, "rowData"))

setReplaceMethod("rowRanges", c("GenomicFiles", "GenomicRangesORGRangesList"),
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
### dim, dimnames
###

setMethod(dim, "GenomicFiles",
    function(x)
{
    c(length(rowRanges(x)), nrow(colData(x)))
})

setMethod(dimnames, "GenomicFiles",
    function(x)
{
    list(names(rowRanges(x)), rownames(colData(x)))
})

setReplaceMethod("dimnames", c("GenomicFiles", "list"),
    function(x, value)
{
    if (length(value) != 2)
        stop("'value' must be a list of length 2")
    names(rowRanges(x)) <- value[[1]]
    rownames(colData(x)) <- value[[2]]
    x
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
        initialize(x, rowData=rowRanges(x)[i,])
)

setMethod("[", c("GenomicFiles", "missing", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (is.character(j))
        j <- match(j, colnames(x))
    if (any(is.na(j)))
        stop("subscript 'j' out of bounds")
    initialize(x, 
               colData=colData(x)[j,,drop=FALSE],
               files=files(x)[j])
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
        files=files(x)[j],
        rowData=rowRanges(x)[i,],
        colData=colData(x)[j,,drop=FALSE])
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Show
###

setMethod(show, "GenomicFiles", 
    function(object) 
{
    cat(class(object), "object with", 
        paste(dim(object), c("ranges", "files:"), collapse=" and "),
        "\n")
    cat("files:", paste(BiocGenerics:::selectSome(basename(files(object))), 
        collapse=", "), "\n")
    cat("detail: use files(), rowRanges(), colData(), ...",
        "\n")
})
