

setMethod(writePhredQual, "FastqQuality", function(object, file, mode="w", ...) {
    if (length(file) != 1)
        sprintf("UserArgumentMismatch:'%s' must be '%s'",
                       "file", "character(1)")
    if (file.exists(file) && mode != "a")
        sprintf("UserArgumentMismatch:file '%s' exists, but mode is not 'a'",
                       file)
    file <- path.expand(file)
    ## FIXME: different quality types
    max_width <- max(c(unique(width(names(sread(object)))),
                       unique(width(quality(object)))))
    .Call("write_phred_quality", id(object), 
          quality(object), file, mode, max_width)
    invisible(length(object))
})
