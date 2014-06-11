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
    signature="ranges")

setMethod(reduceByRange, "GRangesList", .reduceByRange) 

setMethod(reduceByRange, "GRanges", 
    function(ranges, files, MAPPER, REDUCER, iterate=FALSE, ..., init)
        .reduceByRange(as(ranges, "List"), files, MAPPER, 
                       REDUCER, iterate, ..., init=init)
) 
