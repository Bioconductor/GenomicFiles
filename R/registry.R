### =========================================================================
### Utilities for creating and searching a 'file type' registry 
### =========================================================================

.fileTypeRegistry <- new.env(parent=emptyenv())

registerFileType <- function(type, package, regex)
{
    .fileTypeRegistry[[regex]] <- list(package=package, type=type)
    invisible(.fileTypeRegistry[[regex]])
}

findTypeRegistry <- function(fnames)
{
    regexes <- ls(.fileTypeRegistry)
    for (regex in regexes)
        if (all(grepl(regex, fnames)))
            return(regex)
    stop("unknown file type ", paste(sQuote(fnames), collapse=", "))
}

makeFileType <- function(fnames, ..., regex=findTypeRegistry(fnames))
{
    nmspc <- getNamespace(.fileTypeRegistry[[regex]]$package)
    type <- .fileTypeRegistry[[regex]]$type
    FUN <- get(type, nmspc)
    do.call(FUN, list(fnames,  ...))
}
