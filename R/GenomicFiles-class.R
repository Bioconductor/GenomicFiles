### =========================================================================
### GenomicFiles class
### =========================================================================

setGeneric(".validity", function(object) standardGeneric(".validity"))

setClass("GenomicFiles",
    contains="RangedSummarizedExperiment",
    representation(
        files="ANY"
    ),
    prototype(
        files=character()),
    validity=.validity
)


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
    function(rowRanges, files, ...) standardGeneric("GenomicFiles"),
    signature=c("rowRanges", "files"))

### Combine the new parallel slots with those of the parent class.
setMethod("parallelSlotNames", "GenomicFiles",
    function(x) c("files", callNextMethod())
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

setMethod(GenomicFiles, c("GenomicRangesORGRangesList", "character"),
   function(rowRanges, files, colData=DataFrame(), metadata=list(), ...)
{
    if (length(files)) {
        if (is.null(nms <- names(files))) {
            nms <- basename(files)
            names(files) <- nms
        }
        if (missing(colData))
            colData <- DataFrame(row.names=nms)
        else
            rownames(colData) <- nms 
    }
    new("GenomicFiles", rowRanges=rowRanges, 
                        files=files, 
                        colData=colData, 
                        metadata=metadata, ...)
})

setMethod(GenomicFiles, c("GenomicRangesORGRangesList", "List"),
   function(rowRanges, files, colData=DataFrame(), metadata=list(), ...)
{
    if (length(files)) {
        if (is.null(nms <- names(files)))
            stop("'List' of files must be named")

        if (missing(colData))
            colData <- DataFrame(row.names=basename(nms))
        else
            rownames(colData) <- basename(nms)
    }
    new("GenomicFiles", rowRanges=rowRanges, 
                        files=files, 
                        colData=colData, 
                        metadata=metadata, ...)
})

setMethod(GenomicFiles, c("GenomicRangesORGRangesList", "list"),
    function(rowRanges, files, ...)
{
    GenomicFiles(rowRanges, as(files, "List"), ...)
})

setMethod(GenomicFiles, c("missing", "ANY"),
    function(rowRanges, files, ...)
{
    GenomicFiles(GRanges(), files, ...)
})

setMethod(GenomicFiles, c("missing", "missing"),
    function(rowRanges, files, ...)
{
    GenomicFiles(GRanges(), character(), ...)
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
        nms <- basename(value)
    colData <- colData(x)
    rownames(colData) <- nms
    initialize(x, colData=colData, files=value)
})

setReplaceMethod("files", c("GenomicFiles", "List"),
    function(x, ..., value)
{
    if (is.null(nms <- names(value)))
        nms <- value
    colData <- colData(x)
    rownames(colData) <- nms
    initialize(x, colData=colData, files=value)
})

setReplaceMethod("colData", c("GenomicFiles", "DataFrame"),
    function(x, ..., value)
{
    if (length(files(x)) != nrow(value))
        stop("'length(files(x))' must equal 'nrow(value)'")
    files <- files(x)
    names(files) <- rownames(value)
    initialize(x, colData=value, files=files)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting
###

setMethod("[", c("GenomicFiles", "ANY", "missing"),
    function(x, i, j, ..., drop=TRUE)
        initialize(x, rowRanges=rowRanges(x)[i,])
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
        rowRanges=rowRanges(x)[i,],
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
