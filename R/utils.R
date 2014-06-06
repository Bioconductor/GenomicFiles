
.reduce <-
    function(..., BY=c("range", "file"))
{
    BY <- match.arg(BY)
    if (BY == "range")
        reduceByRange(...)
    else
        reduceByFile(...)
}

## TODO: promote to generic / methods
.summarizeView <-
    function(X, MAPPER, REDUCER, ..., BY=c("range", "file"))
{
    BY <- match.arg(BY)
    assay <- .reduce(fileRange(X), fileList(X), 
                     MAPPER, REDUCER, ..., BY=BY)
    if (BY == "range") {
        ## reformat list structure
        if (is(assay, "list")) {
            assay <- bplapply(seq_along(fileList(X)), 
                function(i) sapply(assay, "[", i))
            assay <- simplify2array(assay)
        } else {
            assay <- t(assay)
        }
    } else {
        assay <- simplify2array(assay)
    }
    SummarizedExperiment(assay, rowData=fileRange(X),
                         colData=DataFrame(filePath=path(fileList(X))))
}
