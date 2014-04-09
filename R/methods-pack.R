### =========================================================================
### pack and unpack
### =========================================================================

.pack <- function(x, pack_type, max_len, max_inter_range_len, 
                  width, density_ratio, pack_map=FALSE)
{
    original <- x
    x_order <- order(seqnames(x))
    x <- splitAsList(x, as.character(seqnames(x)))
    x_range <- range(x)
    i_range <- psetdiff(x_range, x)
    i_dist <- i_range[width(i_range) > max_inter_range_len]
    range_LT_max_len <- unlist(width(x_range), use.names=FALSE) < max_len
    is_distant <- elementLengths(i_dist) > 0
    is_dense <- any(sum(width(reduce(x))) > density_ratio * width(x_range))
    is_long <- any(width(x) > max_len)

    LENGTH <- is_long & pack_type %in% c("length", "all")
    DENSITY <- is_dense & !is_distant & !is_long &
               range_LT_max_len & pack_type %in% c("density", "all")
    DISTANCE <- is_distant & !is_long &
                pack_type %in% c("distance", "all")

    if (!any(LENGTH | DENSITY | DISTANCE)) {
        if (pack_map)
            return(new("Hits"))
        else
            return(original)
    }
    ## long -> tile
    if (any(LENGTH)) {
        tile_fun <- function(i, width) {
            require(IRanges) ## for windows
            unlist(tile(i, width=width), use.names=FALSE)
        }
        x[LENGTH] <- GRangesList(bplapply(x[LENGTH], tile_fun, width)) 
    }
    ## dense -> max range
    if (any(DENSITY))
        x[DENSITY] <- x_range[DENSITY]
    ## distant -> create bins
    if (any(DISTANCE))
        x[DISTANCE] <- psetdiff(x_range[DISTANCE], i_dist[DISTANCE]) 

    if (pack_map)
        findOverlaps(original, unlist(x, use.names=FALSE))
    else
        unlist(x, use.names=FALSE)
}

setMethod("pack", "GRanges",
    function(x, ..., pack_type = c("all", "density", "distance", "length"),
             extraction_type = c("NumericList", "RleList"),
             max_len = 1e9, max_inter_range_len = 1e7, width = 1e7, 
             density_ratio = 0.25) 
{ 
    pack_type <- match.arg(pack_type)
    extraction_type <- match.arg(extraction_type)
    if (extraction_type == "RleList")
        return(x)
    if (width > max_len)
        stop("'width' must be less than 'max_len'")
    if (density_ratio < 0 | density_ratio > 1)
        stop("'density_ratio' must be >= 0 and <= 1")

    .pack(x, pack_type=pack_type, max_len=max_len, 
          max_inter_range_len=max_inter_range_len, width=width,
          density_ratio=density_ratio, ...)
})

setMethod("unpack", c("GRanges", "NumericList"),
    function(x, y, ...) 
{
    if (is.null(metadata(y)$ranges))
        stop("metadata(y)$ranges must contain corresponding IRanges")

    map <- pack(x, pack_map=TRUE, ...)
    s_hits <- subjectHits(map)
    q_hits <- queryHits(map)
    y_range <- metadata(y)$ranges
    x_group <- split(s_hits, q_hits)
    x_group_min_start <- sapply(x_group, 
        function(i, y_range) min(start(y_range[i])), y_range)
    x_shift <- shift(x, - (x_group_min_start - 1))

    ## FIXME: use unlist/relist to combine tiled y ranges
    y_group <- unname(split(s_hits, q_hits))
    y_order <- vector("list", length(x_shift))
    if (any(idx <- elementLengths(y_group) > 1)) {
        FUN <- function(i, y) {
            require(IRanges) ## for windows
            NumericList(unlist(y[i]))
        }
        y_order <- do.call(c, bplapply(y_group, FUN, y))
        x_shift <- ranges(splitAsList(x_shift, seq_along(y_order)))
    } else {
        y_order <- y[unique(s_hits)]
        x_shift <- ranges(splitAsList(x_shift, match(s_hits, s_hits)))
    }
    res <- do.call(c, bpmapply(extractList, y_order, x_shift, USE.NAMES=FALSE))
    names(res) <- seqnames(x)
    res
})

setMethod("unpack", c("GRanges", "RleList"),
    function(x, y, ...) 
{
    ## FIXME: explore findRange() and findRun()
    y[x]
})
