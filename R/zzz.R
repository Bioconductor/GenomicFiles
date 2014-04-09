.onLoad <- function(libname, pkgname)
{
    registerFileType("FaFileList", "Rsamtools", "\\.fa$")
    registerFileType("FaFileList", "Rsamtools", "\\.fasta$")
    registerFileType("BamFileList", "Rsamtools", "\\.bam$")
    registerFileType("BigWigFileList", "rtracklayer", "\\.bw$")
}
