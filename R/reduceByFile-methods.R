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
    if (missing(REDUCER) && iterate)
        iterate <- FALSE

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
                result <- REDUCER(list(result, mapped), ...)
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
    function(ranges, files, MAPPER, REDUCER, iterate=TRUE, ..., init)
        standardGeneric("reduceByFile"),
    signature=c("ranges", "files")
)

setMethod(reduceByFile, c("GRangesList", "ANY"), 
    function(ranges, files, MAPPER, REDUCER, iterate=TRUE,
             summarize=FALSE, ..., init) {
        lst <- .reduceByFile(ranges, files, MAPPER, REDUCER, iterate, 
                             ..., init)
        if (summarize && missing(REDUCER))
            SummarizedExperiment(simplify2array(lst), rowData=ranges,
                                 colData=DataFrame(filePath=files))
        else
            lst
    }
) 

setMethod(reduceByFile, c("GRanges", "ANY"), 
    function(ranges, files, MAPPER, REDUCER, iterate=TRUE, 
             summarize=FALSE, ..., init) {
        lst <- .reduceByFile(as(ranges, "List"), files, MAPPER, REDUCER,
                             iterate, ..., init=init)
        if (summarize && missing(REDUCER))
            SummarizedExperiment(simplify2array(lst), rowData=ranges,
                                 colData=DataFrame(filePath=files))
        else
            lst
    }
) 

setMethod(reduceByFile, c("GenomicFiles", "missing"), 
    function(ranges, files, MAPPER, REDUCER, iterate=TRUE, 
             summarize=FALSE, ..., init) {
        reduceByFile(rowData(ranges)[[1]], GenomicFiles::files(ranges),
                     MAPPER, REDUCER, iterate, summarize, ..., init=init)
    }
)
