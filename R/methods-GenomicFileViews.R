### =========================================================================
### GenomicFileViews (VIRTUAL) 
### =========================================================================

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
### reduceByFile and reduceByRange 
###
.reduce <-
    function(..., BY=c("range", "file"))
{
    BY <- match.arg(BY)
    if (BY == "range")
        reduceByRange(...)
    else
        reduceByFile(...)
}

setMethod(reduceByRange, "GenomicFileViews",
    function(X, MAP, REDUCE, ..., init, ITERATE=FALSE)
{
    NO_REDUCE <- missing(REDUCE)
    if(ITERATE && NO_REDUCE) 
       stop("'ITERATE' must be FALSE when 'REDUCE' is missing")
    if (!NO_REDUCE) 
        REDUCE <- match.fun(REDUCE)
    MAP <- match.fun(MAP)
    initmiss <- missing(init)
    grl <- unname(splitAsList(fileRange(X), seq_along(fileRange(X))))

    bplapply(grl, function(grange, ...) {
        if (ITERATE) {
            result <- if (initmiss) {
                MAP(fileList(X)[[1]], grange, ...)
            } else init
            for (i in seq_along(fileList(X))[-1]) {
                mapped <- MAP(fileList(X)[[i]], grange, ...)
                result <- REDUCE(result, mapped, ..., RANGE=grange)
            }
            result
        } else {
            mapped <- lapply(fileList(X), MAP, RANGE=grange, ...)
            if (NO_REDUCE)
                mapped
            else
                REDUCE(mapped, ..., RANGE=grange)
        }
    }, ...)
})

setMethod(reduceByFile, "GenomicFileViews",
    function(X, MAP, REDUCE, ..., init, ITERATE=FALSE)
{
    NO_REDUCE <- missing(REDUCE)
    if(ITERATE && NO_REDUCE) 
       stop("'ITERATE' must be FALSE when 'REDUCE' is missing")
    if (!NO_REDUCE) 
        REDUCE <- match.fun(REDUCE)
    MAP <- match.fun(MAP)
    initmiss <- missing(init)
    grl <- unname(splitAsList(fileRange(X), seq_along(fileRange(X))))

    bplapply(fileList(X), function(fl, ...) {
        if (ITERATE) {
            result <- if (initmiss) {
                MAP(fl, grl[[1]], ...)
            } else init
            for (i in seq_along(grl)[-1]) {
                mapped <- MAP(fl, grl[[i]], ...)
                result <- REDUCE(result, mapped, ..., FILE=fl)
            }
            result
        } else {
            mapped <- lapply(grl, MAP, FILE=fl, ...)
            if (NO_REDUCE)
                mapped
            else
                REDUCE(mapped, ..., FILE=fl)
        }
    }, ...)
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
