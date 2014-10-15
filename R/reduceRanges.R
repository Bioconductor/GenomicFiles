### =========================================================================
### reduceRanges
### =========================================================================

reduceRanges <- function(ranges, files, MAP, REDUCE, ..., iterate=TRUE, init) {
    if (is(ranges, "GenomicFiles")) {
        files <- GenomicFiles::files(ranges)
        ranges <- rowData(ranges)
    }
    if (!is(ranges, "GRanges") && !is(ranges, "GRangesList"))
        stop("'ranges' must be GRanges or GRangesList")

    reduceByRange(ranges, list(files), MAP, REDUCE, ..., 
                  summarize=FALSE, iterate=iterate)
}
