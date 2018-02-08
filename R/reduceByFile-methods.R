### =========================================================================
### reduceByFile
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generic and methods
###

.reduceByFile <- function(ranges, files, MAP, REDUCE, ..., iterate, init)
{
    if (!is(files, "character") && !is(files, "List"))
        stop("'files' must be character vector or List of filenames")
    if (missing(REDUCE) && iterate)
        iterate <- FALSE
    if (missing(REDUCE))
        REDUCE <- NULL
    if (missing(init))
        init <- NULL

    ## files sent to workers
    bplapply(files, function(file, ranges, MAP, REDUCE, ..., iterate, init) {
        if (iterate) {
            result <- 
                if (is.null(init))
                    MAP(ranges[[1]], file, ...)
                else init
            for (i in seq_along(ranges)[-1]) {
                mapped <- MAP(ranges[[i]], file, ...)
                result <- REDUCE(list(result, mapped), ...)
            }
            result
        } else {
            mapped <- lapply(ranges, MAP, file, ...)
            if (is.null(REDUCE))
                mapped
            else
                REDUCE(mapped, ...)
        }
    }, ..., ranges=ranges, MAP=MAP, REDUCE=REDUCE, iterate=iterate, init=init)
}

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
        if (summarize && !missing(REDUCE))
            warning("'summarize' set to FALSE when REDUCE is provided")
        if (summarize && missing(REDUCE))
            SummarizedExperiment(SimpleList(list(data=simplify2array(lst))), 
                                 rowRanges=ranges, 
                                 colData=DataFrame(filePath=files))
        else
            lst
    }
) 

setMethod(reduceByFile, c("GRanges", "ANY"), 
    function(ranges, files, MAP, REDUCE, ..., summarize=FALSE,
             iterate=TRUE, init) {
        lst <- .reduceByFile(as(ranges, "CompressedGRangesList"), files, MAP,
                             REDUCE, ..., iterate=iterate, init=init)
        if (summarize && !missing(REDUCE))
            warning("'summarize' set to FALSE when REDUCE is provided")
        if (summarize && missing(REDUCE))
            SummarizedExperiment(SimpleList(list(data=simplify2array(lst))), 
                                 rowRanges=ranges,
                                 colData=DataFrame(filePath=files))
        else
            lst
    }
) 

setMethod(reduceByFile, c("GenomicFiles", "missing"), 
    function(ranges, files, MAP, REDUCE, ..., summarize=FALSE,
             iterate=TRUE, init) {
        reduceByFile(rowRanges(ranges), GenomicFiles::files(ranges),
                     MAP, REDUCE, ..., summarize=summarize,
                     iterate=iterate, init=init)
    }
)
