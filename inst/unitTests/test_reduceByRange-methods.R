fl <- system.file("extdata", "ex1.bam", package="Rsamtools") 
gr <- GRanges(c("seq1", "seq2", "seq2"), IRanges(1:3, width=100))
gf <- GenomicFiles(gr, c(one=fl, two=fl))

MAP = function(RANGE, FILE, ..., param=Rsamtools::ScanBamParam()) {
    Rsamtools::bamWhich(param) <- RANGE 
    Rsamtools::countBam(FILE, param=param)
}

## reduceByRange
test_reduceByRange_MAP <- function()
{
    ans <- reduceByRange(gf, MAP=MAP)
    checkIdentical(length(ans), 3L)
    checkIdentical(unname(lengths(ans)), c(2L, 2L, 2L))

    ans <- reduceByRange(gr, c(one=fl, two=fl), MAP=MAP)
    checkIdentical(length(ans), 3L)
    checkIdentical(unname(lengths(ans)), c(2L, 2L, 2L))

    ## summarize = TRUE
    ans <- reduceByRange(gf, MAP=MAP, summarize=TRUE)
    checkTrue(is(ans, "SummarizedExperiment"))
    checkIdentical(names(assays(ans)), 'data')
    checkIdentical(dim(assays(ans)$data), c(3L, 2L))
}

test_reduceByRange_MAP_REDUCE <- function()
{
    REDUCE = function(MAPPED, ...) {
        head(MAPPED, 1)
    }

    ans <- reduceByRange(gf, MAP=MAP, REDUCE=REDUCE)
    checkIdentical(length(ans), 3L)
    checkIdentical(unname(elementNROWS(ans)), c(1L, 1L, 1L))
}

## reduceRanges

## TBD
