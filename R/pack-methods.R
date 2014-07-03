### =========================================================================
### pack methods 
### =========================================================================
 
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generic and methods
###

setGeneric("pack", function(x, ...)
    standardGeneric("pack"),
    signature="x")

setMethod("pack", "GRanges",
    function(x, ..., range_len=1e9, inter_range_len=1e7)
        .pack(x, range_len=range_len, inter_range_len=inter_range_len)
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Helpers
###

isPacked <- function(x, ...)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    if (is(x@partitioning, "PartitioningMap"))
        TRUE 
    else
        FALSE 
}

.pack <- function(x, range_len, inter_range_len) 
{
    if (length(x) == 0)
        return(x) 

    ## order
    o <- order(x)
    as.character(seqnames(x))
    if (is.unsorted(o))
        x_grl <- splitAsList(x[o], seqnames(x)[o])
    else
        x_grl <- splitAsList(x, seqnames(x))

    ## identify 'long' and 'distant'
    long <- which(width(unlist(x_grl, use.names=FALSE)) > range_len)
    long_minus1 <- long - 1L
    long_minus1 <- long_minus1[long_minus1 > 0L]
    irange <- unname(psetdiff(range(x_grl), x_grl))
    irange_max <- irange[width(irange) > inter_range_len]
    irange_idx <- elementLengths(irange_max) > 0
    distant <- integer()
    if (any(irange_idx))
        distant <- sapply(irange_max[irange_idx], 
            function(i, xx) 
                follow(i, xx, ignore.strand=TRUE),
            xx=unlist(x_grl, use.names=FALSE))

    ends <- c(distant, long_minus1, long, end(PartitioningByEnd(x_grl)))
    x_grl@partitioning <- PartitioningMap(x=sort(unique(ends)), order(o)) 
    x_grl
}

