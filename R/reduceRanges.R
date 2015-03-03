### =========================================================================
### reduceRanges
### =========================================================================

reduceRanges <- function(ranges, files, MAP, REDUCE, ..., init) {
    if (is(ranges, "GenomicFiles")) {
        files <- GenomicFiles::files(ranges)
        ranges <- rowRanges(ranges)
    }
    if (!is(ranges, "GRanges") && !is(ranges, "GRangesList"))
        stop("'ranges' must be GRanges or GRangesList")

    .reduceByRange(ranges, list(files), MAP, REDUCE, ..., 
                   iterate=FALSE)
}
