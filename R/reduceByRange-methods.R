### =========================================================================
### Queries across files (reduceByRange) 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generic and methods
###

.reduceByRange <- function(ranges, files, MAP, REDUCE, ..., iterate, init)
{
    if (!is(files, "character") && !is(files, "List"))
        stop("'files' must be character vector or List of filenames")
    if (missing(REDUCE) && iterate)
        iterate <- FALSE
    if (missing(REDUCE))
        REDUCE <- NULL
    if (missing(init))
        init <- NULL

    ## ranges sent to workers
    bplapply(ranges, function(elt, files, MAP, REDUCE, ..., iterate, init) {
        require(GenomicRanges)
        if (iterate) {
            result <- if (is.null(init)) {
                MAP(elt, files[[1]], ...)
            } else init
            for (i in seq_along(files)[-1]) {
                mapped <- MAP(elt, files[[i]], ...)
                result <- REDUCE(list(result, mapped), ...)
            }
            result
        } else {
            mapped <- lapply(files, function(f) MAP(elt, f, ...))
            if (is.null(REDUCE))
                mapped
            else
                REDUCE(mapped, ...)
        }
    }, ..., files=files, MAP=MAP, REDUCE=REDUCE, iterate=iterate, init=init)
}

setGeneric("reduceByRange", 
    function(ranges, files, MAP, REDUCE, ..., iterate=TRUE, init)
        standardGeneric("reduceByRange"),
    signature=c("ranges", "files")
)

setMethod(reduceByRange, c("GRangesList", "ANY"), 
    function(ranges, files, MAP, REDUCE, ..., summarize=FALSE,
             iterate=TRUE, init) {
        lst <- .reduceByRange(ranges, files, MAP, REDUCE, ...,
                              iterate=iterate)
        if (summarize && !missing(REDUCE))
            warning("'summarize' set to FALSE when REDUCE is provided")
        if (summarize && missing(REDUCE)) {
            lst <- bplapply(seq_along(files), 
                function(i) sapply(lst, "[", i))
            SummarizedExperiment(SimpleList(list(data=simplify2array(lst))), 
                                 rowData=ranges,
                                 colData=DataFrame(filePath=files))
        } else {
            lst
        }
    }
)

setMethod(reduceByRange, c("GRanges", "ANY"), 
    function(ranges, files, MAP, REDUCE, ..., summarize=FALSE,
             iterate=TRUE, init) {
        lst <- .reduceByRange(as(ranges, "List"), files, MAP, REDUCE,
                              ..., iterate=iterate) 
        if (summarize && !missing(REDUCE))
            warning("'summarize' set to FALSE when REDUCE is provided")
        if (summarize && missing(REDUCE)) {
            lst <- bplapply(seq_along(files), 
                function(i) sapply(lst, "[", i))
            SummarizedExperiment(SimpleList(list(data=simplify2array(lst))), 
                                 rowData=ranges,
                                 colData=DataFrame(filePath=files))
        } else {
            lst
        }
    }
) 

setMethod(reduceByRange, c("GenomicFiles", "missing"), 
    function(ranges, files, MAP, REDUCE, ..., summarize=FALSE,
             iterate=TRUE, init) {
        reduceByRange(rowData(ranges), GenomicFiles::files(ranges),
                      MAP, REDUCE, ..., iterate=iterate)
    }
)
