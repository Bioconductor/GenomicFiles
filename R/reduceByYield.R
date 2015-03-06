### =========================================================================
### reduceByYield (iterate through files by chunk)
### =========================================================================

.reduceByYield_iterate <-
    function(X, YIELD, MAP, REDUCE, DONE, ..., parallel, init)
{
    if (parallel) {
        ITER <- function() {
            if(DONE(value <- YIELD(X, ...)))
                NULL
            else
               value
        }
        result <- bpiterate(ITER, FUN=MAP, REDUCE=REDUCE, ...)
    } else {
        result <- if (missing(init)) {
            data <- YIELD(X, ...)
            if (DONE(data))
                return(list())
            MAP(data, ...)
        } else
            init

        repeat {
            if(DONE(data <- YIELD(X, ...)))
                break
            value <- MAP(data, ...)
            result <- REDUCE(result, value)
        }
    }
    result
}
 
.reduceByYield_all <-
    function(X, YIELD, MAP, REDUCE, DONE, ..., parallel)
{
    if (parallel) {
        ITER <- function() {
            if(DONE(value <- YIELD(X, ...)))
                NULL
            else
               value
        }
        result <- bpiterate(ITER, FUN=MAP, ...)
    } else {
        result <- bpiterate(ITER, FUN=MAP, ..., BPPARAM=SerialParam())
    }
    REDUCE(result)
}

## REDUCE and init are never NULL; init can be missing
reduceByYield <-
    function(X, YIELD, MAP = identity, REDUCE = `+`, 
             DONE = function(x) is.null(x) || length(x) == 0L, 
             ..., parallel=FALSE, iterate=TRUE, init)
{
    if  (missing(REDUCE)) 
        REDUCE <- if (iterate) c else identity
    if (!iterate && !missing(init))
        warning("'init' ignored when iterate == FALSE")

    if (!isOpen(X)) {
        open(X)
        on.exit(close(X))
    }
    if (iterate)
        .reduceByYield_iterate(X, YIELD, MAP, REDUCE, DONE, ...,
                               parallel=parallel, init=init)
    else
        .reduceByYield_all(X, YIELD, MAP, REDUCE, DONE,
                           ..., parallel=parallel)
}

REDUCEsampler <-
    function(sampleSize=1000000, verbose=FALSE)
{
    tot <- 0L
    function(x, y, ...) {
        if (length(x) < sampleSize)
            stop("expected yield of at least sampleSize=", sampleSize)

        if (tot == 0L) {
            ## first time through
            tot <<- length(x)
            x <- x[sample(length(x), sampleSize)]
        }
        yld_n <- length(y)
        tot <<- tot + yld_n

        if (verbose)
            message("REDUCEsampler total=", tot)

        keep <- rbinom(1L, min(sampleSize, yld_n), yld_n / tot)
        i <- sample(sampleSize, keep)
        j <- sample(yld_n, keep)
        x[i] <- y[j]
        x
    }
}
