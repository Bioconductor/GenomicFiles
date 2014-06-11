### =========================================================================
### Perform queries within files (reduceByFile)
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers 
###

.reduceByFile <- function(ranges, files, MAPPER, REDUCER, iterate, ..., init)
{
    if (!is(files, "character") && !is(files, "List"))
        stop("'files' must be character vector or List of filenames")
    if (missing(REDUCER) && iterate==TRUE)
        stop("when 'iterate=TRUE' 'REDUCER' must be specified")

    NO_REDUCER <- missing(REDUCER)
    ## files sent to workers
    bplapply(files, function(file, ..., init) {
        if (iterate) {
            result <- 
                if (missing(init))
                    MAPPER(ranges[[1]], file, ...)
                else init
            for (i in seq_along(ranges)[-1]) {
                mapped <- MAPPER(ranges[[i]], file, ...)
                result <- REDUCER(result, mapped, ...)
            }
            result
        } else {
            mapped <- lapply(ranges, MAPPER, file, ...)
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

setGeneric("reduceByFile", 
    function(ranges, files, MAPPER, REDUCER, iterate=FALSE, ..., init)
    standardGeneric("reduceByFile"),
    signature="ranges")

setMethod(reduceByFile, "GRangesList", .reduceByFile) 

setMethod(reduceByFile, "GRanges", 
    function(ranges, files, MAPPER, REDUCER, iterate=FALSE, ..., init)
        .reduceByFile(as(ranges, "List"), files, MAPPER, 
                      REDUCER, iterate, ..., init=init)
) 
