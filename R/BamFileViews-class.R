### =========================================================================
### BamFileViews methods
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity
###

setMethod(.validity, "BamFileViews", 
    function(object) {
        msg <- NULL
        if (!is(object@fileList, "BamFileList"))
            msg <- "'fileList' must be a 'BamFileList'"
        if (is.null(msg)) TRUE else msg
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructors 
###

setMethod(BamFileViews, "missing", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(path(fileList)))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    stop("'fileList' must be character() or BamFileList")
})

setMethod(BamFileViews, "BamFileList", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(path(fileList)))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    new("BamFileViews", ..., 
        fileList=fileList,
        fileSample=fileSample, fileRange=fileRange,
        fileExperiment=fileExperiment, yieldSize=yieldSize, 
        .views_on_file=.views_on_file)
})

setMethod(BamFileViews, "character", 
          function(fileList,
                   fileSample=DataFrame(row.names=
                     make.unique(basename(fileList))),
                   fileRange=GRanges(),
                   fileExperiment=list(),
                   yieldSize=NA_integer_,
                   .views_on_file=new.env(parent=emptyenv()), ...)
{
    new("BamFileViews", ..., 
        fileList=BamFileList(lapply(fileList, BamFile, ...)),
        fileSample=fileSample, fileRange=fileRange,
        fileExperiment=fileExperiment, yieldSize=yieldSize, 
        .views_on_file=.views_on_file)
})

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### summarizeOverlaps(), countBam(), scanBam() 
###

.compare_ranges <- function(file, param)
{ 
    ranges1 <- fileRange(file)
    ranges2 <- bamWhich(param)
    if (length(ranges1))  {
      if (length(ranges2) && !identical(ranges1, ranges2))
          warning("'fileRange(reads)' and 'bamWhich(param)' differ; ",
                  "using fileRange(reads)")
      bamWhich(param) <- ranges1
    }
    param
}

.so_BFV <- 
    function(features, reads, mode, ignore.strand=FALSE, ..., 
             inter.feature=TRUE, singleEnd=TRUE, fragments=FALSE,
             param=ScanBamParam())
{
    param <- .compare_ranges(reads, param)
    summarizeOverlaps(features, fileList(reads), mode,
        ignore.strand=ignore.strand, inter.feature=inter.feature, 
        singleEnd=singleEnd, fragments=fragments, param=param, ...)
}
setMethod("summarizeOverlaps", c("GRanges", "BamFileViews"), .so_BFV)
setMethod("summarizeOverlaps", c("GRangesList", "BamFileViews"), .so_BFV)


.cb_BFV <- function(file, index=file, ..., param=ScanBamParam())
{
    param <- .compare_ranges(file, param)
    bplapply(fileList(file), countBam, index=fileList(index), ..., param=param)
}
setMethod(countBam, "BamFileViews", .cb_BFV)


.sc_BFV <- function(file, index=file, ..., 
                    param=ScanBamParam(what=scanBamWhat()))
{
    param <- .compare_ranges(file, param)
    bplapply(fileList(file), scanBam, index=fileList(index), ..., param=param)
}
setMethod(scanBam, "BamFileViews", .sc_BFV)

