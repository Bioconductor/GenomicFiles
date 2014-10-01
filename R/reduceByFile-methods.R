### =========================================================================
### Queries within files (reduceByFile)
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers 
###

.reduceByFile <- function(ranges, files, MAP, REDUCE, ..., iterate, init)
{
    if (!is(files, "character") && !is(files, "List"))
        stop("'files' must be character vector or List of filenames")
    if (missing(REDUCE) && iterate)
        iterate <- FALSE

    NO_REDUCE <- missing(REDUCE)
    ## files sent to workers
    bplapply(files, function(file, ..., init) {
        if (iterate) {
            result <- 
                if (missing(init))
                    MAP(ranges[[1]], file, ...)
                else init
            for (i in seq_along(ranges)[-1]) {
                mapped <- MAP(ranges[[i]], file, ...)
                result <- REDUCE(list(result, mapped), ...)
            }
            result
        } else {
            mapped <- lapply(ranges, MAP, file, ...)
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

setGeneric("reduceByFile", 
    function(ranges, files, MAP, REDUCE, ..., iterate=TRUE, init)
        standardGeneric("reduceByFile"),
    signature=c("ranges", "files")
)

setMethod(reduceByFile, c("GRangesList", "ANY"), 
    function(ranges, files, MAP, REDUCE, ..., summarize=FALSE,
             iterate=TRUE, init) {
        lst <- .reduceByFile(ranges, files, MAP, REDUCE, 
                             ..., iterate=iterate, init=init)
        if (summarize && missing(REDUCE))
            SummarizedExperiment(simplify2array(lst), rowData=ranges,
                                 colData=DataFrame(filePath=files))
        else
            lst
    }
) 

setMethod(reduceByFile, c("GRanges", "ANY"), 
    function(ranges, files, MAP, REDUCE, ..., summarize=FALSE,
             iterate=TRUE, init) {
        lst <- .reduceByFile(as(ranges, "List"), files, MAP, 
                             REDUCE, ..., iterate=iterate, init=init)
        if (summarize && missing(REDUCE))
            SummarizedExperiment(simplify2array(lst), rowData=ranges,
                                 colData=DataFrame(filePath=files))
        else
            lst
    }
) 

setMethod(reduceByFile, c("GenomicFiles", "missing"), 
    function(ranges, files, MAP, REDUCE, ..., summarize=FALSE,
             iterate=TRUE, init) {
        reduceByFile(rowData(ranges), GenomicFiles::files(ranges),
                     MAP, REDUCE, ..., summarize=summarize,
                     iterate=iterate, init=init)
    }
)
