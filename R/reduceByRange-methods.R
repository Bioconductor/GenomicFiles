### =========================================================================
### Perform queries across files (reduceByRange) 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers 
###

.reduceByRange <- function(ranges, files, MAPPER, REDUCER, iterate, ..., init)
{
    if (!is(files, "character") && !is(files, "List"))
        stop("'files' must be character vector or List of filenames")
    if (missing(REDUCER) && iterate==TRUE)
        stop("when 'iterate=TRUE' 'REDUCER' must be specified")

    NO_REDUCER <- missing(REDUCER)
    ## ranges sent to workers
    bplapply(ranges, function(elt, ..., init) {
        if (iterate) {
            result <- if (missing(init)) {
                MAPPER(elt, files[[1]], ...)
            } else init
            for (i in seq_along(files)[-1]) {
                mapped <- MAPPER(elt, files[[i]], ...)
                result <- REDUCER(result, mapped, ...)
            }
            result
        } else {
            mapped <- lapply(files, function(f) MAPPER(elt, f, ...))
            if (NO_REDUCER)
                mapped
            else
                REDUCER(mapped, ...)
        }
    }, ...)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generic and methods
###

setGeneric("reduceByRange", 
    function(ranges, files, MAPPER, REDUCER, iterate=FALSE, ..., init)
        standardGeneric("reduceByRange"),
    signature=c("ranges", "files")
)

setMethod(reduceByRange, c("GRangesList", "ANY"), 
    function(ranges, files, MAPPER, REDUCER, iterate=FALSE, 
             summarize=FALSE, ..., init) {
        if (summarize && !missing(REDUCER))
            stop("'summarize' must be FALSE when REDUCER is provided")
        lst <- .reduceByRange(ranges, files, MAPPER, REDUCER, iterate,
                              ..., init=init)
        if(summarize) {
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
    function(ranges, files, MAPPER, REDUCER, iterate=FALSE, 
             summarize=FALSE, ..., init) {
        if (summarize && !missing(REDUCER))
            stop("'summarize' must be FALSE when REDUCER is provided")
        lst <- .reduceByRange(as(ranges, "List"), files, MAPPER,
                              REDUCER, iterate, ..., init=init) 
        if(summarize) {
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
    function(ranges, files, MAPPER, REDUCER, iterate=FALSE, 
             summarize=FALSE, ..., init) {
        if (summarize && !missing(REDUCER))
            stop("'summarize' must be FALSE when REDUCER is provided")
        reduceByRange(rowData(ranges)[[1]], GenomicFiles::files(ranges),
                      MAPPER, REDUCER, iterate, summarize, ..., init=init)
    }
)
