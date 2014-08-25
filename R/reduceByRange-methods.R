### =========================================================================
### Queries across files (reduceByRange) 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers 
###

.reduceByRange <- function(ranges, files, MAP, REDUCE, ..., iterate, init)
{
    if (!is(files, "character") && !is(files, "List"))
        stop("'files' must be character vector or List of filenames")
    if (missing(REDUCE) && iterate)
        iterate <- FALSE

    NO_REDUCE <- missing(REDUCE)
    ## ranges sent to workers
    bplapply(ranges, function(elt, ..., init) {
        if (iterate) {
            result <- if (missing(init)) {
                MAP(elt, files[[1]], ...)
            } else init
            for (i in seq_along(files)[-1]) {
                mapped <- MAP(elt, files[[i]], ...)
                result <- REDUCE(list(result, mapped), ...)
            }
            result
        } else {
            mapped <- lapply(files, function(f) MAP(elt, f, ...))
            if (NO_REDUCE)
                mapped
            else
                REDUCE(mapped, ...)
        }
    }, ...)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generic and methods
###

setGeneric("reduceByRange", 
    function(ranges, files, MAP, REDUCE, ..., iterate=TRUE, init)
        standardGeneric("reduceByRange"),
    signature=c("ranges", "files")
)

setMethod(reduceByRange, c("GRangesList", "ANY"), 
    function(ranges, files, MAP, REDUCE, ..., summarize=FALSE,
             iterate=TRUE, init) {
        lst <- .reduceByRange(ranges, files, MAP, REDUCE, ...,
                              iterate=iterate, init=init)
        if (summarize && missing(REDUCE)) {
            lst <- bplapply(seq_along(files), 
                function(i) sapply(lst, "[", i))
            SummarizedExperiment(simplify2array(lst), rowData=ranges,
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
                              ..., iterate=iterate, init=init) 
        if (summarize && missing(REDUCE)) {
            lst <- bplapply(seq_along(files), 
                function(i) sapply(lst, "[", i))
            SummarizedExperiment(simplify2array(lst), rowData=ranges,
                                 colData=DataFrame(filePath=files))
        } else {
            lst
        }
    }
) 

setMethod(reduceByRange, c("GenomicFiles", "missing"), 
    function(ranges, files, MAP, REDUCE, ..., summarize=FALSE,
             iterate=TRUE, init) {
        reduceByRange(rowData(ranges)[[1]], GenomicFiles::files(ranges),
                      MAP, REDUCE, ..., summarize=summarize,
                      iterate=iterate, init=init)
    }
)
