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
    ans <- reduceFiles(gf, MAP=MAP)
    checkIdentical(length(ans), 2L)
    checkIdentical(unname(elementLengths(ans)), c(1L, 1L))

    ans <- reduceFiles(gr, c(one=fl, two=fl), MAP=MAP)
    checkIdentical(length(ans), 2L)
    checkIdentical(unname(elementLengths(ans)), c(1L, 1L))

    ans <- reduceFiles(gr, c(one=fl, two=fl), MAP=MAP, nchunk=2)
    checkIdentical(length(ans), 2L)
    checkIdentical(unname(elementLengths(ans)), c(2L, 2L))
}

 test_reduceFiles_MAP_REDUCE <- function()
{
    REDUCE <- function(mapped) do.call(rbind, mapped)
    ans <- reduceFiles(gf, MAP=MAP, REDUCE=REDUCE, nchunk=2)
    checkIdentical(length(ans), 2L)
    checkIdentical(unname(elementLengths(ans)), c(1L, 1L))
}

