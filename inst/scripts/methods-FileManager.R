### =========================================================================
### FileManager objects
### =========================================================================

### Tasks to manage:
### - (split) divide problem for parallel eval
### - (apply / MAP) function executed on workers
### - (combine / REDUCE) collate results from workers 

### Thoughts:
### - Parallel exection:
###   Should be optional (large data, shared memory).
###
### - chunk/yield/iterate:
###   These are are independent of parallel and appropriate 
###   to include in apply/MAP step.
###
### - I/O considerations:
###   Potential problems with parallel reduceByRange() and 
###   single ff or hdf5 file (all hit same file, same time).
###
### - run byFile but reduce byRange:
###   Extraction of multiple ranges is more efficient
###   using a param (reduceByFile approach). If memory isn't 
###   a problem, we may consider a fast method to combine
###   elements of the list returned by bplapply() across
###   list elements. This would allow the user to run byFile
###   but reduce byRange.


###  Reference Class

.FM <- setRefClass("FileManager",
    fields = list(
        files = "ANY",
        chunkSplit = "function",
        chunkApply = "function",
        chunkCombine = "function",
        parallel = "BiocParallelParam"),
    methods = list(
        initialize = function(
            parallel = bpparam(), yieldSize=integer(), 
            chunkSplit = function(x=.self$files) {
                if (is(x, "List")) x
                else split(x, seq_along(x))
            },
            chunkCombine = function(x) unlist(x, use.names=FALSE), ...) {
                callSuper(...)
                .self$parallel=parallel
                .self$chunkSplit=chunkSplit 
                .self$chunkCombine=chunkCombine
                .self
        },
        run = function(...) {
            res <- bplapply(.self$chunkSplit(), .self$chunkApply, ...)
            .self$chunkCombine(res)
        }))


## case 1: whole file 
fl <- system.file("extdata", "ex1.bam", package="Rsamtools")
fm1 <- .FM$new(files=c(fl, fl, fl), chunkApply=countBam)

## case 2: ranges
param <- ScanBamParam(which=GRanges("seq1", IRanges(1, 20)), what="qname") 
fm2 <- .FM$new(files=c(fl, fl, fl), chunkApply=scanBam) 
fm2$run(param=param)

## case 3: BamFile single yield
bfl <- BamFileList(c(fl, fl, fl), yieldSize=50L) 
fm3 <- .FM$new(files=bfl, chunkApply=scanBam, chunkCombine=return)
fm3$run(param=ScanBamParam(what="qname"))

## case 4: roll your own yield
fun <- function(i, ...) {
    open(i)
    on.exit(close(i))
    ans <- NULL
    while (length(res <- scanBam(i, ...)[[1]]$qname))
        ans <- c(ans, length(res))
    ans
}
fm4 <- .FM$new(files=bfl, chunkApply=fun)
fm4$run()


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Standard S4 
###

setClass("FileManager",
    representation(
        files="character",
        chunkSplit="function",
        chunkApply="function",
        chunkCombine="function",
        parallel="BiocParallelParam"),
    prototype(
        parallel=bpparam()),
    validity=.validity)

## validity would confirm required methods

## constructor
FileManager <- 
    function(files, chunkApply, 
             chunkSplit = function(object) length(files),
             chunkExtract = function(object, i) files[[i]],
             chunkCombine = function(x) unlist(x, use.names=FALSE),
             parallel=bpparam(), ...)
{
    new("FileManager", files, chunkApply, chunkSplit, chunkExtract, 
        chunkCombine, parallel)
}

