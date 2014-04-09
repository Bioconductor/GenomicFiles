### =========================================================================
### .FileList 
### =========================================================================

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Validity. 
###

setMethod(.validity, ".FileList", 
    function(object) {
        msg <- NULL
        if (!c("path", "index") %in% names(object))
            msg <- c(msg, 
                     "names(fileList(object)) must include 'path' and 'index'")
        if (is.null(msg)) TRUE else msg
    }
)
