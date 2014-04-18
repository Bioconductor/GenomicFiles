### =========================================================================
### pack and unpack methods
### =========================================================================
 
.pack <- function(x, range_len, inter_range_len) 
{
    if (length(x) == 0)
        return(x) 

    ## order
    o <- order(x)
    if (is.unsorted(o))
        x_grl <- split(x[o], as.character(seqnames(x)[o]))
    else
        x_grl <- split(x, as.character(seqnames(x))) 

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

setMethod("pack", "GRanges",
    function(x, ..., range_len=1e9, inter_range_len=1e7)
        .pack(x, range_len=range_len, inter_range_len=inter_range_len)
)

isPacked <- function(x, ...)
{
    if (!is(x, "GRangesList"))
        stop("'x' must be a GRangesList object")
    if (is(x@partitioning, "PartitioningMap"))
        TRUE 
    else
        FALSE 
}

## handle results from *lapply()
setMethod("unpack", c("list", "GRangesList"),
    function(flesh, skeleton, ...)
        unpack(List(flesh), skeleton, ...)
) 

setMethod("unpack", c("List", "GRangesList"),
    function(flesh, skeleton, ...) 
{
    if (!isPacked(skeleton))
        stop("'flesh' must be a packed object")

    mo <- mapOrder(skeleton@partitioning)
    if (is(flesh, "RleList"))
        do.call(c, flesh)[mo]
    else
        unlist(flesh, use.names=FALSE)[mo]
})
