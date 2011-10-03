
setMethod(.sffValidity, "SffReads", function(object) {
    msg <- NULL
    lenq <- length(quality(object))
    lens <- length(sread(object))
    if (lenq != lens) {
        txt <- sprintf("sread and quality length mismatch: %d %d",
                       lenq, lens)
        msg <- c(msg, txt)
    }
    if (!all(width(quality(object)) == width(sread(object)))) {
        txt <- sprintf("some sread and quality widths differ")
        msg <- c(msg, txt)
    }
    if (is.null(msg)) TRUE else msg
})

setMethod(SffReads, c("DNAStringSet", "QualityScore", "IRanges","IRanges","character"),
          function(sread, quality, qualityClip, adapterClip, mode, ...)
{
    new("SffReads", sread=sread, quality=quality, 
        qualityClip=qualityClip, adapterClip=adapterClip, clipMode=mode, ...)
})

setMethod(SffReads, c("DNAStringSet", "QualityScore", "IRanges","IRanges","missing"),
          function(sread, quality, qualityClip, adapterClip, mode, ...)
{
    new("SffReads", sread=sread, quality=quality,
        qualityClip=qualityClip, adapterClip=adapterClip, clipMode="Full", ...)
})

setMethod(SffReads, c("DNAStringSet", "QualityScore", "IRanges","missing","missing"),
          function(sread, quality, qualityClip, adapterClip, mode, ...)
{
	emptyIR <- IRanges(start=integer(length(sread)),width=integer(length(sread)))
    new("SffReads", sread=sread, quality=quality,
        qualityClip=qualityClip, adapterClip=emptyIR, clipMode="Full", ...)
})

setMethod(SffReads, c("DNAStringSet", "QualityScore", "missing","IRanges","missing"),
          function(sread, quality, qualityClip, adapterClip, mode, ...)
{
	emptyIR <- IRanges(start=integer(length(sread)),width=integer(length(sread)))
    new("SffReads", sread=sread, quality=quality,
        qualityClip=emptyIR, adapterClip=adapterClip, clipMode="Full", ...)
})

setMethod(SffReads, c("DNAStringSet", "QualityScore", "missing","missing","missing"),
          function(sread, quality, qualityClip, adapterClip, mode, ...)
{
	emptyIR <- IRanges(start=integer(length(sread)),width=integer(length(sread)))
    new("SffReads", sread=sread, quality=quality,
        qualityClip=emptyIR, adapterClip=emptyIR, clipMode="Full", ...)
})

setMethod(SffReads, c("DNAStringSet", "BStringSet", "IRanges", "IRanges", "character"),
    function(sread, quality, qualityClip, adapterClip, mode, ...,
             qualityType="FastqQuality")
{
        quality <- FastqQuality(quality)
	    new("SffReads", sread=sread, quality=quality,
	        qualityClip=qualityClip, adapterClip=adapterClip, clipMode=mode, ...)
})

setMethod(SffReads, c("DNAStringSet", "BStringSet", "IRanges", "IRanges", "missing"),
    function(sread, quality, qualityClip, adapterClip, mode, ...,
             qualityType="FastqQuality")
{
        quality <- FastqQuality(quality)
	    new("SffReads", sread=sread, quality=quality,
	        qualityClip=emptyIR, adapterClip=emptyIR, clipMode="Full", ...)
})

setMethod(SffReads, c("missing", "missing", "missing","missing","missing"),
          function(sread, quality, qualityClip, adapterClip, mode, ...)
{
    new("SffReads")
})

#ShortRead:::.make_getter("sread")

#setMethod(sread, "SffReads",
#         function(object, ...) slot(object, "sread"))

##setMethod(quality, "SffReads",
##          function(object, ...) slot(object, "quality"))

setMethod(id, "SffReads", function(x) names(sread(x)))

setMethod(length, "SffReads", function(x) length(sread(x)))

setMethod(width, "SffReads", function(x) width(sread(x)))


setMethod(writeFastq, "SffReads", function(object, file, mode="w", ...) {
    if (length(file) != 1)
        .throw(SRError("UserArgumentMismatch", "'%s' must be '%s'",
                       "file", "character(1)"))
    if (file.exists(file) && mode != "a")
        .throw(SRError("UserArgumentMismatch",
                       "file '%s' exists, but mode is not 'a'",
                       file))
    file <- path.expand(file)
    ## FIXME: different quality types
    max_width <- max(c(unique(width(names(sread(object)))),
                       unique(width(sread(object))),
                       unique(width(quality(object)))))
    .Call(.write_fastq, names(sread(object)), sread(object),
          quality(quality(object)), file, mode, max_width)
    invisible(length(object))
})

## show
setMethod(show, "SffReads", function(object) {
    callNextMethod()
    wd <- sort(unique(width(object)))
    if (length(wd)>2) wd <- paste(range(wd), collapse="..")
    cat("length:", length(object), "reads; width:", wd, "basepair\n")
})
## detail
setMethod(detail, "SffReads", function(x, ...) {
    cat("class: ", class(x), "\n")
    cat("\nsread:\n")
    show(sread(x))
})
