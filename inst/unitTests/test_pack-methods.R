FUN_IL <- function(i) IntegerList(as.list(start(i)))
FUN_RL <- function(i) RleList(as.list(start(i)))

.unpack <- function(pk, gr)
{
    dat <- lapply(pk, FUN_IL)
    upk <- unpack(dat, pk)
    checkTrue(length(gr) == length(upk))
    checkIdentical(start(gr), unlist(upk))

    dat <- lapply(pk, FUN_RL)
    upk <- unpack(dat, pk)
    checkTrue(length(gr) == length(upk))
    checkIdentical(start(gr), runValue(unlist(upk)))
}

test_pack_unpack_no_op <- function()
{
    ## empty
    gr <- GRanges()
    checkIdentical(pack(gr), GRanges())
    grl <- GRangesList()
    grl@partitioning <- PartitioningMap()
    checkIdentical(unpack(IntegerList(), grl), integer())

    ## no re-order or re-group
    gr <- GRanges("chr1", IRanges(1:5*5, width=3))
    pk <- pack(gr)
    checkIdentical(ranges(unlist(pk)), ranges(gr))
    .unpack(pk, gr)
}

test_pack_unpack_order <- function()
{
    gr <- GRanges(c(rep("chr4", 3), "chr1", "chr1"), 
        IRanges(c(10, 1, 100, 5, 2), width=1))
    pk <- pack(gr)
    checkTrue(length(pk) == 2L)
    .unpack(pk, gr)
}

test_pack_unpack_distant <- function()
{
    gr1 <- GRanges("chr1", IRanges(c(1, 5, 30000, 30005), width=3))
    pk1 <- pack(gr1, inter_range_len=1000)
    checkTrue(length(pk1) == 2L)
    pm <- pk1@partitioning
    checkIdentical(width(pm), c(2L, 2L))
    checkIdentical(mapOrder(pm), as.integer(1:4))

    gr2 <- GRanges("chr1", IRanges(c(1, 30000, 30005), width=3))
    pk2 <- pack(gr2, inter_range_len=1000)
    checkTrue(length(pk2) == 2L)
    pm <- pk2@partitioning
    checkIdentical(width(pm), c(1L, 2L))
    checkIdentical(mapOrder(pm), as.integer(1:3))

    gr3 <- GRanges("chr1", IRanges(c(1, 5, 30000), width=3))
    pk3 <- pack(gr3, inter_range_len=1000)
    checkTrue(length(pk3) == 2L)
    pm <- pk3@partitioning
    checkIdentical(width(pm), c(2L, 1L))
    checkIdentical(mapOrder(pm), as.integer(1:3))

    .unpack(pk1, gr1)
    .unpack(pk2, gr2)
    .unpack(pk3, gr3)
}

test_pack_unpack_long <- function()
{
    gr1 <- GRanges("chr2", IRanges(c(20, 10, 1), width=c(5, 2000, 5)))
    pk1 <- pack(gr1, range_len=1000)
    checkTrue(length(pk1) == 3L)
    pm <- pk1@partitioning
    checkIdentical(width(pm), c(1L, 1L, 1L))
    checkIdentical(mapOrder(pm), c(3L, 2L, 1L))

    gr2 <- GRanges("chr2", IRanges(c(1, 10, 20), width=c(2000, 5, 5)))
    pk2 <- pack(gr2, range_len=1000)
    checkTrue(length(pk2) == 2L)
    pm <- pk2@partitioning
    checkIdentical(width(pm), c(1L, 2L))
    checkIdentical(mapOrder(pm), as.integer(1:3))

    gr3 <- GRanges("chr2", IRanges(c(1, 10, 20), width=c(5, 5, 2000)))
    pk3 <- pack(gr3, range_len=1000)
    checkTrue(length(pk3) == 2L)
    pm <- pk3@partitioning
    checkIdentical(width(pm), c(2L, 1L))
    checkIdentical(mapOrder(pm), as.integer(1:3))

    .unpack(pk1, gr1)
    .unpack(pk2, gr2)
    .unpack(pk3, gr3)
}
