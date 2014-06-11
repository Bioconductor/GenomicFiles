### =========================================================================
### pack and unpack methods
### =========================================================================
 
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Generic and methods
###

setGeneric("unpack", function(flesh, skeleton, ...)
    standardGeneric("unpack"),
    signature=c("flesh", "skeleton"))

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
