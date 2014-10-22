### =========================================================================
### reduceFiles
### =========================================================================

reduceFiles <- function(ranges, files, MAP, REDUCE, ..., init) {
    if (is(ranges, "GenomicFiles")) {
        files <- GenomicFiles::files(ranges)
        ranges <- rowData(ranges)
    }
    if (!is(ranges, "GRanges") && !is(ranges, "GRangesList"))
        stop("'ranges' must be GRanges or GRangesList")

    .reduceByFile(list(ranges), files, MAP, REDUCE, ...,
                 summarize=FALSE, iterate=FALSE)
}
