### pack

test_pack_no_pack <- function()
{
    gr <- GRanges("chr1", IRanges(1:5*5, width=3))
    ans1 <- pack(gr, pack_type="distance")
    checkIdentical(ans1, gr)
    ans2 <- pack(gr, pack_type="distance", pack_map=TRUE)
    checkIdentical(ans2, new("Hits"))

    ans1 <- pack(gr, pack_type="length")
    checkIdentical(ans1, gr)
    ans2 <- pack(gr, pack_type="length", pack_map=TRUE)
    checkIdentical(ans2, new("Hits"))

    gr <- GRanges("chr1", IRanges(c(1, 50), width=3))
    ans1 <- pack(gr, pack_type="density")
    checkIdentical(ans1, gr)
    ans2 <- pack(gr, pack_type="density", pack_map=TRUE)
    checkIdentical(ans2, new("Hits"))
}

test_pack_density <- function()
{
    gr <- GRanges("chr1", IRanges(1:5*5, width=3))

    ans0 <- range(gr)
    ans1 <- pack(gr)
    ans2 <- pack(gr, pack_type="density")
    checkIdentical(ans0, ans1)
    checkIdentical(ans1, ans2)

    ans0 <- findOverlaps(gr, pack(gr))
    ans1 <- pack(gr, pack_map=TRUE)
    checkIdentical(ans0, ans1)
}

test_pack_distance <- function()
{
    gr1 <- GRanges("chr4", 
        IRanges(c(1, 10, 20, 30, 40, 75), width=c(rep(10, 5), 25)))
    gr2 <- GRanges("chr3", IRanges(c(1, 5, 30000, 30005), width=3))
    gr <- suppressWarnings(c(gr1, gr2))

    ans1 <- pack(gr)
    ans2 <- pack(gr, pack_type="distance")
    checkIdentical(ans1, ans2)

    ans0 <- findOverlaps(gr, pack(gr))
    ans1 <- pack(gr, pack_map=TRUE)
    checkIdentical(ans0, ans1)
}

test_pack_length <- function()
{
    gr <- GRanges("chr2", IRanges(c(1, 10), width=c(5, 2000)))
    max_len <- 1000

    ans1 <- pack(gr, max_len=max_len, width=500)
    ans2 <- pack(gr, pack_type="length", max_len=max_len, width=500)
    checkIdentical(ans1, ans2)

    ans0 <- findOverlaps(gr, pack(gr, max_len=max_len, width=500))
    ans1 <- pack(gr, max_len=max_len, width=500, pack_map=TRUE)
    checkIdentical(ans0, ans1)

    ans0 <- findOverlaps(gr, pack(gr, max_len=1000, width=250))
    ans3 <- pack(gr, max_len=1000, width=250, pack_map=TRUE)
    checkIdentical(ans0, ans3)
    checkTrue(!identical(ans1, ans3))
}

### unpack

fl <- system.file("tests", "test.bw", package = "rtracklayer")
test_unpack_density <- function()
{
    x <- GRanges("chr2", IRanges(1:10*5, width=4)) 
    y <- import(fl, which=pack(x), as="NumericList")
    ans1 <- unpack(x, y)
    checkTrue(length(x) == length(ans1))
} 

test_unpack_distance <- function()
{
    x <- GRanges(c(rep("chr2", 3), rep("chr19", 3)), 
        IRanges(c(1, 5, 450, 1450, 1750, 1755), width=3)) 
    y <- import(fl, which=pack(x), as="NumericList")
    ans1 <- unpack(x, y)
    checkTrue(length(x) == length(ans1))
    checkTrue(sum(width(x)) == sum(elementLengths(ans1)))
}

test_unpack_length <- function()
{
    max_len <- width <- 500
    x <- GRanges(c("chr2", "chr2", "chr19"), 
        IRanges(c(1, 105, 1), width=c(100, 2000, 10))) 
    p <- pack(x, max_len=max_len, width=width)
    checkIdentical((x %in% p), c(TRUE, FALSE, TRUE))

    y <- import(fl, which=p, as="NumericList")
    ans1 <- unpack(x, y, max_len=max_len, width=width)
    checkTrue(length(x) == length(ans1))
    checkTrue(sum(width(x)) == sum(elementLengths(ans1)))
}

test_unpack_distance_and_length <- function()
{
    ## length takes priority
    max_len <- width <- 500
    max_inter_range <- 200
    x <- GRanges(c("chr2", "chr2"), 
        IRanges(c(1, 500), width=c(10, 1000))) 
    p <- pack(x, max_len=max_len, width=width)
    checkIdentical((x %in% p), c(TRUE, FALSE))

    y <- import(fl, which=p, as="NumericList")
    ans1 <- unpack(x, y, max_len=max_len, width=width)
    checkTrue(length(x) == length(ans1))
    checkTrue(sum(width(x)) == sum(elementLengths(ans1)))
}








