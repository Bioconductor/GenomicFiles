.summarizeView <-
    function(X, MAP, REDUCE, ..., BY=c("range", "file"))
{
    BY <- match.arg(BY)
    assay <- .reduce(X, MAP, REDUCE, ..., BY=BY)
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
