extdata <- system.file(package="GenomicFiles", "extdata")
files <- dir(extdata, pattern="^CEUtrio.*bgz$", full=TRUE)
names(files) <- sub(".*_([0-9XY]+).*", "\\1", basename(files))
seqinfo <- as(readRDS(file.path(extdata, "seqinfo.rds")), "Seqinfo")
smps <- samples(VariantAnnotation::scanVcfHeader(files[1]))
colData <- DataFrame(row.names=smps)

test_VcfStack_construction <- function() {
    ## empty constructor
    checkTrue(validObject(VcfStack()))

    ## named files
    checkException(validObject(VcfStack(unname(files))))

    ## all files must exist
    checkException(VcfStack(tempfile()))

    ## constructor with files only
    checkTrue(validObject(VcfStack(files)))

    ## constructor with files and seqinfo
    checkTrue(validObject(VcfStack(files, seqinfo)))

    ## constructor with files and wrong seqinfo
    nm <- names(files)
    names(files)[1] <- "BreakThis"
    checkException(VcfStack(files, seqinfo))
    names(files) <- nm

    ## constructor with files, seqinfo, and colData
    checkTrue(validObject(VcfStack(files, seqinfo, colData)))

    ## constructor with files and colData
    checkTrue(validObject(VcfStack(files, colData=colData)))

    ## constructor with files and wrong colData
    checkException(VcfStack(
        files,
        colData=DataFrame(row.names=c("Break", "This", "Now"))))

    ## check override check - bad!
    checkTrue(validObject(VcfStack(files,
           colData=DataFrame(row.names=c("Break","This", "Now")),
                                   check =FALSE)))

    ## constructor with seqinfo and colData
    checkTrue(validObject(VcfStack(seqinfo=seqinfo, colData=colData)))

    ## constructor with seqinto
    checkTrue(validObject(VcfStack(seqinfo=seqinfo)))

    ## constructor with colData
    checkTrue(validObject(VcfStack(colData=colData)))
}

test_RangedVcfStack_construction <- function() {

    ## empty constructor
    checkTrue(validObject(RangedVcfStack()))

    ## constructor with VcfStack object
    checkTrue(validObject(RangedVcfStack(VcfStack(files))))

    ## constructor with valid rowRanges object
    checkTrue(validObject(RangedVcfStack(
        VcfStack(files),
        rowRanges=GRanges(c("7:1-100000000","X:1-100000000")))))

    ## constructor with invalid rowRanges object
    checkException(RangedVcfStack(
        VcfStack(files),
        rowRanges=GRanges(
            c("7:1-100000000", "X:1-100000000", "19:1-100000000"))))
}

test_RangedVcfStack_construction_2 <- function() {
    stack <- VcfStack(files)
    gr0 <- GRanges(seqinfo(stack))[rownames(stack)]

    checkTrue(validObject(RangedVcfStack(stack, gr0)))

    gr1 <- gr0[1:2]
    checkTrue(validObject(RangedVcfStack(stack, gr1)))

    gr2 <- GRanges(seqinfo(stack))[c("1", "2")]
    checkException(RangedVcfStack(stack, gr2))
}

test_RangedVcfStack_seqinfo <- function() {
    rstack <- RangedVcfStack(VcfStack(files))
    value0 <- seqinfo(rstack)

    ## valid update of seqinfo -- reduce seqlevels
    value <- value0[rownames(rstack)]
    seqinfo(rstack) <- value
    checkIdentical(seqinfo(rstack), value)
    checkIdentical(seqinfo(rowRanges(rstack)), value)

    ## fail to drop seqlevels currently in use
    value <- value0[setdiff(seqnames(value0), rownames(rstack))]
    checkException({
        seqinfo(rstack) <- value
    })
}

test_VcfStack_subsetting <- function() {

    ## default data is 7 files x 3 samples object
    stack <- VcfStack(files, seqinfo)

    ## test empty subsetting
    checkTrue(all(dim(stack[])==c(7,3)))

    ## test numeric subsetting
    checkTrue(all(dim(stack[1:3,])==c(3,3)))
    checkTrue(all(dim(stack[,2])==c(7,1)))
    checkTrue(all(dim(stack[1:3,2])==c(3,1)))

    ## test character subsetting
    ## default object names
    ## files: "11" "20" "21" "22" "7"  "X"  "Y"
    ## samples: "NA12878" "NA12891" "NA12892"
    checkTrue(all(dim(stack[c("X", "11"),])==c(2,3)))
    checkTrue(all(dim(stack[,c("NA12878","NA12891")])==c(7,2)))
    checkTrue(all(dim(stack[c("X", "11"),"NA12892"])==c(2,1)))

    ## test mix numeric and character subsetting
    checkTrue(all(dim(stack[c("X", "11"),1:2])==c(2,2)))
    checkTrue(all(dim(stack[1:3,"NA12892"])==c(3,1)))

    ## test GRange object subsetting
    checkTrue(all(dim(stack[GRanges("20:862167-62858306")])==c(1,3)))
    checkTrue(all(dim(stack[GRanges("20:862167-62858306"),1])==c(1,1)))
    checkTrue(all(dim(stack[GRanges("20:862167-62858306"),"NA12891"])==c(1,1)))

    ## errors if out of bounds or value not found
    checkException(stack[4:8,])
    checkException(stack[,2:4])
    checkException(stack[c("X", "19")],)
    checkException(stack[,c("NA12878","NOTFOUND")])
    checkException(stack[GRanges("19:1-235466666")])

}


test_RangedVcfStack_subsetting <- function(){

    # VcfStack object with 7 files and 3 samples
    # GRanges object with 2 ranges and 0 metadata columns
    Rstack <- RangedVcfStack(VcfStack(files, seqinfo),
                             GRanges(c("7:1-159138000", "X:1-155270560")))

    # empty subset
    checkIdentical(dim(Rstack[,]), dim(Rstack))

    # check sample subsetting
    checkIdentical(dim(Rstack[, 1]), c(7L, 1L))
    checkIdentical(dim(Rstack[,c(TRUE, FALSE, TRUE)]), c(7L, 2L))
    checkIdentical(dim(Rstack[,"NA12891"]), c(7L, 1L))

    # check file subsetting and updating GRanges object
    checkIdentical(dim(Rstack[1,]), c(1L, 3L))
    checkIdentical(length(seqnames(rowRanges(Rstack[1,]))), 0L)
    checkIdentical(dim(Rstack["7",]), c(1L, 3L))
    checkIdentical(as.character(seqnames(rowRanges(Rstack["7",]))), "7")
    gr <- GRanges(c("X:1-100000"))
    checkIdentical(as.character(seqnames(rowRanges(Rstack[gr,]))), "X")
    checkException(Rstack[GRanges(c("X:1-100000", "13:1-100000")),])
}

test_VcfStack_replaceFiles <- function(){

    stack <- VcfStack(files)
    checkTrue(class(files(stack)) == "VcfFileList")
    checkTrue(all(dim(stack) == c(7L, 3L)))

    # replace with character
    files(stack) = files[1]
    checkTrue(all(dim(stack) == c(1L, 3L)))

    # replace with VcfFileList
    files(stack) = VariantAnnotation::VcfFileList(files[1:3])
    checkTrue(all(dim(stack) == c(3L, 3L)))
}

test_VcfStack_readVcfStack <- function(){

    stack <- VcfStack(files)

    # all files
    temp = readVcfStack(stack)
    checkTrue(all(dim(temp) == c(1000L, 3L)))

    # test read by numeric and read by character
    temp1 = readVcfStack(stack, 1)
    temp2 = readVcfStack(stack, names(files(stack[1])))
    checkTrue(all(dim(temp1) == dim(temp2)))

    # test read by GRange
    gr = GRanges("11:1-100")
    temp3 = readVcfStack(stack, gr)
    checkTrue(all(dim(temp3) == c(0L, 3L)))
    gr = GRanges(paste0("11:1-", seqlengths(seqinfo(stack))[levels(seqnames(gr))]))
    temp3 = readVcfStack(stack, gr)
    checkTrue(all(dim(temp1) == dim(temp3)))

    # test read by range out of bounds
    checkException(readVcfStack(stack, GRanges("11:1-1000000000")))

    # check multiple
    temp4 = readVcfStack(stack, 3)
    temp5 = readVcfStack(stack,c(1,3))
    checkTrue((dim(temp1)[1] + dim(temp4)[1]) == dim(temp5)[1])


    # test RangedVcfStack
    gr = GRanges("11:1-135006516")
    Rstack = RangedVcfStack(stack, gr)
    temp6 = readVcfStack(Rstack)
    checkTrue(all(dim(temp6) == dim(temp3)))

    gr = GRanges(c("11:1-135006516", "21:1-48129895"))
    Rstack = RangedVcfStack(stack, gr)
    temp7 = readVcfStack(Rstack)
    checkTrue(all(dim(temp7) == dim(temp5)))
}
