## Used with GenomicFileViews class 

.reduce <-
    function(..., BY=c("range", "file"))
{
    BY <- match.arg(BY)
    if (BY == "range")
        reduceByRange(...)
    else
        reduceByFile(...)
}

.summarizeView <-
    function(X, MAPPER, ..., BY=c("range", "file"))
{
    BY <- match.arg(BY)
    assay <- .reduce(fileRange(X), fileList(X), 
                     MAPPER, ..., BY=BY)
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
    SummarizedExperiment(assay, rowRanges=fileRange(X),
                         colData=DataFrame(filePath=path(fileList(X))))
}
