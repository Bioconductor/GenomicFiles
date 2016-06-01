extdata <- system.file(package="GenomicFiles", "extdata")
files <- dir(extdata, pattern="^CEUtrio.*bgz$", full=TRUE)
files2 <- files
names(files) <- sub(".*_([0-9XY]+).*", "\\1", basename(files))
seqinfo <- as(readRDS(file.path(extdata, "seqinfo.rds")), "Seqinfo")
smps <- samples(VariantAnnotation::scanVcfHeader(files[1]))
colData <- DataFrame(row.names=smps)

test_VcfStack_construction <- function() {
    ## empty constructor
    checkTrue(validObject(VcfStack()))

    ## named files
    checkException(validObject(VcfStack(files2)))

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
    checkException(VcfStack(files, colData=DataFrame(row.names=c("Break", "This", "Now"))))

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
    checkTrue(validObject(RangedVcfStack(VcfStack(files), rowRanges=GRanges(c("7:1-100000000","X:1-100000000")))))

    ## constructor with invalid rowRanges object 
    checkException(RangedVcfStack(VcfStack(files), rowRanges=GRanges(c("7:1-100000000","X:1-100000000", "19:1-100000000"))))

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
    Rstack <- RangedVcfStack(VcfStack(files, seqinfo),GRanges(c("7:1-159138000","X:1-155270560")))

    # empty subset 
    checkTrue(all(dim(Rstack[])==c(2,3)))
   
    # check sample subsetting 
    checkTrue(all(dim(Rstack[,1])==c(2,1)))
    checkTrue(all(dim(Rstack[,c(T,F,T)])==c(2,2)))
    checkTrue(all(dim(Rstack[,"NA12891"])==c(2,1)))

    # check file subsetting - always use rowRanges
    checkTrue(all(dim(Rstack[1,])==c(2,3)))
    checkTrue(all(dim(Rstack["11",])==c(2,3)))

}

