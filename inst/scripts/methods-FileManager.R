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

## use yield_reduce 
## iterators pkg icount 
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

## lazy split
## load_balancing
## should have more tasks than nodes but not more than
## by a factor of 100 (or even 10?)

## sequential, non-scheduled, no knowledge of length

    } else {
        ## pre-allocate result list to retain order
        sx <- seq_along(X)
        res <- vector("list", length(sx))
        names(res) <- names(X)

        ## start as many jobs as there are cores 
        ## don't need to wait for all jobs to finish
        ## can access elements of 'jobs' individually
        jobid <- seq_len(cores)
        jobs <- lapply(jobid,
                       function(i) mcparallel(FUN(X[[i]], ...),
                                              mc.set.seed = mc.set.seed,
                                              silent = mc.silent))
        jobsp <- processID(jobs)
        has.errors <- 0L
        complete <- 0L 
        current <- lazySplit(length(X)) ## iterator
        curr <- nextElem(current)
        ## advance curr to number of cores
        while (curr < cores) curr <- nextElem(current)
        while (length(i <- nextElem(current))) {
            s <- selectChildren(jobs, 0.5)
            if (is.null(s)) break   # no children -> no hope
            if (is.integer(s)) { # one or more children finished
                for (ch in s) {
                    ji <- which(jobsp == ch)[1]
                    ci <- jobid[ji]
                    r <- readChild(ch)
                    if (is.raw(r)) {
                        child.res <- unserialize(r)
                        if (inherits(child.res, "try-error"))
                            has.errors <- has.errors + 1L
                        ## a NULL assignment would remove it from the list
                        if (!is.null(child.res)) res[[ci]] <- child.res
                        complete <- complete + 1
                    } else {
                        complete <- complete + 1
                        if (length(curr <- nextElem(current))) {
                            # spawn a new job
                            nexti <- curr 
                            jobid[ji] <- nexti
                            jobs[[ji]] <- mcparallel(FUN(X[[nexti]], ...),
                                                     mc.set.seed=mc.set.seed,
                                                     silent=mc.silent)
                            jobsp[ji] <- processID(jobs[[ji]])
                        }
                    }
                }
            }
        }






