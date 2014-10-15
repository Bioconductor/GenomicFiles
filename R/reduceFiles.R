### =========================================================================
### reduceFiles
### =========================================================================

.chunkIndex <- function(rows, nchunk, ...) 
{
    if (missing("nchunk")) {
        if (rows > 1e8) 
            n <- ceiling(rows / 3)
        else if (rows > 1e6) 
            n <- ceiling(rows / 2)
        else if (rows > 1e5) 
            n <- 1e5 
        else 
            return(list(seq_len(rows)))
    } else {
        if (is.na(nchunk))
            return(list(seq_len(rows)))
        else
            n <- ceiling(rows / nchunk)
    }

    split(seq_len(rows), ceiling(seq_len(rows)/n))
}

reduceFiles <- function(ranges, files, MAP, REDUCE, ..., iterate=TRUE, init) {
    if (is(ranges, "GenomicFiles")) {
        files <- GenomicFiles::files(ranges)
        ranges <- rowData(ranges)
    }
    if (!is(ranges, "GRanges") && !is(ranges, "GRangesList"))
        stop("'ranges' must be GRanges or GRangesList")
    idx <- .chunkIndex(length(ranges), ...)

##    ## Option 1: reduceByRange
##    ## - results must be collated by file before return 
##    ## - only uses 1 core when ranges < 1e5 
#      reduceByRange(relist(ranges, idx), files, MAP, REDUCE, ...,
#                        summarize=FALSE, iterate=iterate)

    ## Option 2: bplapply / bpiterate
    ITER <- function(idx) {
        i <- 0 
        done <- FALSE
  
        function() {
            i <<- i + 1
            if (done)
                return(NULL)
            if (i > length(idx)) {
                done <<- TRUE
                NULL
            } else idx[[i]] 
        }
    }
    MYFUN <- function(file, ITER, idx, MAP, ranges, REDUCE, ...) {
       library(BiocParallel)
       FUN <- function(i, MAP, ranges, file, ...)
           MAP(ranges[i], file, ...)
       res <- bpiterate(ITER(idx), FUN=FUN, MAP, ranges, file, 
                        ..., BPPARAM=SerialParam())
       if (!missing(REDUCE))
           REDUCE(res)
       else
           res
    }
    res <- bplapply(files, FUN=MYFUN, ITER, idx, MAP, ranges, ...)
    if (!missing(REDUCE))
        bplapply(res, REDUCE)
    else
        res
}






