test_GenomicFiles_dimnames <- function() {
    gf <- GenomicFiles(files=c("a", "b"))
    checkIdentical(c("a", "b"), names(files(gf)))
    checkIdentical(list(NULL, c("a", "b")), dimnames(gf))
    colnames(gf) <- c("c", "d")
    checkIdentical(list(NULL, c("c", "d")), dimnames(gf))
    checkIdentical(c("c", "d"), names(files(gf)))
}
