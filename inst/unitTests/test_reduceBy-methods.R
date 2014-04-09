fl <- system.file("extdata", "ex1.bam", package="Rsamtools") 
gr <- GRanges(c("seq1", "seq2", "seq2"), IRanges(1:3, width=100))
bfv <- BamFileViews(c(fl, fl), fileRange=gr)

MAP = function(FILE, RANGE, ..., param) {
    bamWhich(param) <- RANGE 
    countBam(FILE, param=param)
}

REDUCE = function(MAPPED, ...) {
    Reduce("+", MAPPED)
}

## reduceByFile

test_reduceByFile_MAP_only <- function()
{
    ans0 <- reduceByFile(bfv, MAP, param=ScanBamParam())
    checkIdentical(length(ans0), 2L)

    ans1 <- unname(elementLengths(ans0))  
    checkIdentical(ans1, c(3L, 3L))
}

test_reduceByFile_MAP_REDUCE <- function()
{
    ans0 <- reduceByFile(bfv, MAP, REDUCE, param=ScanBamParam())
    checkIdentical(length(ans0), 2L)

    ans1 <- unname(elementLengths(ans0))  
    checkIdentical(ans1, c(1L, 1L))

    ans2 <- unique(unlist(sapply(unname(ans0), "[", "nucleotides")))
    checkIdentical(ans2, 7218)
}

## reduceByRange

test_reduceByRange_MAP_only <- function()
{
    MAP = function(FILE, RANGE, ..., param) {
        bamWhich(param) <- RANGE 
        countBam(FILE, param=param)
    }
    ans0 <- reduceByRange(bfv, MAP, param=ScanBamParam())
    checkIdentical(length(ans0), 3L)

    ans1 <- unname(elementLengths(ans0))  
    checkIdentical(ans1, c(2L, 2L, 2L))
}

test_reduceByRange_MAP_REDUCE <- function()
{
    ans0 <- reduceByRange(bfv, MAP, REDUCE, param=ScanBamParam())
    checkIdentical(length(ans0), 3L)

    ans1 <- unname(elementLengths(ans0))  
    checkIdentical(ans1, c(1L, 1L, 1L))

    ans2 <- unname(unlist(sapply(ans0, "[", "nucleotides")))
    checkIdentical(ans2, c(2744, 5846, 5846))
}
