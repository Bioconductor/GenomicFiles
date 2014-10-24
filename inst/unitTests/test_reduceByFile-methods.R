fl <- system.file("extdata", "ex1.bam", package="Rsamtools") 
gr <- GRanges(c("seq1", "seq2", "seq2"), IRanges(1:3, width=100))
gf <- GenomicFiles(gr, c(one=fl, two=fl))

MAP = function(RANGE, FILE, ..., param=ScanBamParam()) {
    bamWhich(param) <- RANGE 
    countBam(FILE, param=param)
}

REDUCE = function(MAPPED, ...) {
    Reduce("+", MAPPED)
}

## reduceByFile
test_reduceByFile_MAP <- function()
{
    ans <- reduceByFile(gf, MAP=MAP)
    checkIdentical(length(ans), 2L)
    checkIdentical(unname(elementLengths(ans)), c(3L, 3L))

    ans <- reduceByFile(gr, c(one=fl, two=fl), MAP=MAP)
    checkIdentical(length(ans), 2L)
    checkIdentical(unname(elementLengths(ans)), c(3L, 3L))

    ## summarize = TRUE
    ans <- reduceByFile(gf, MAP=MAP, summarize=TRUE)
    checkTrue(is(ans, "SummarizedExperiment"))
    checkIdentical(names(assays(ans)), 'data')
    checkIdentical(dim(assays(ans)$data), c(3L, 2L))
}
 
 test_reduceByFile_MAP_REDUCE <- function()
{
    ans <- reduceByFile(gf, MAP=MAP, REDUCE=REDUCE)
    checkIdentical(length(ans), 2L)
    checkIdentical(unname(elementLengths(ans)), c(1L, 1L))
}

## reduceFiles
test_reduceFiles_MAP <- function()
{
    ## No REDUCE applied, MAP returns 
    ans0 <- reduceFiles(gf, MAP=MAP)
    checkIdentical(length(ans0), 2L)
    checkIdentical(unname(elementLengths(ans0)), c(1L, 1L))
    elts <- lapply(ans0, elementLengths)
    checkIdentical(names(elts), c("one", "two"))
    checkIdentical(unlist(elts, use.names=FALSE), c(3L, 3L))

    ans1 <- reduceFiles(gr, c(one=fl, two=fl), MAP=MAP)
    checkIdentical(ans0, ans1)
    checkIdentical(length(ans1), 2L)
    checkIdentical(unname(elementLengths(ans1)), c(1L, 1L))
}

 test_reduceFiles_MAP_REDUCE <- function()
{
    ## REDUCE applied single time after MAP, simply unlist
    REDUCE <- function(mapped, ...) do.call(rbind, mapped)

    ans <- reduceFiles(gf, MAP=MAP, REDUCE=REDUCE)
    checkIdentical(length(ans), 2L)
    checkIdentical(unname(elementLengths(ans)), c(3L, 3L))
}

