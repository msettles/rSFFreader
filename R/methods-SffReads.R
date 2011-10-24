
## Inspector
setMethod(.sffValidity, "SffReads", function(object) {
#	cat("## SFFreads Validity ###\n")
    msg <- NULL
    lens <- length(object@sread)
    lenqc <- length(object@qualityClip)
    lenac <- length(object@adapterClip)
    if (length(unique(c(lens,lenqc,lenac))) != 1) {
        txt <- sprintf("mismatch length in sread, qualityClip, or adapterClip: %d %d %d",
                       lens,lenqc,lenac)
        msg <- c(msg, txt)
    }
##TODO: need to add check for IRanges out of bounds
	if (!(object@clipMode %in% c("Full","Quality","Raw"))){
		txt <- sprintf("wrong mode type must be one of Full,Quality,Raw")
		msg <- c(msg,txt)
	}
    if (is.null(msg)) TRUE else msg
})

## constructor
"SffReads" <- function(sread, qualityClip, adapterClip,
	clipMode=c("Full", "Quality", "Raw"), header, ...)
{
    if (missing(header)) header = list()
    clipMode = match.arg(clipMode)
    if (missing(qualityClip) | missing(adapterClip))
        emptyIR <- IRanges(start=integer(length(sread)),width=integer(length(sread)))
    if (missing(qualityClip)) qualityClip=emptyIR
    if (missing(adapterClip)) adapterClipj=emtpyIR
    
    new("SffReads", header=header, sread=sread,
        qualityClip=qualityClip, adapterClip=adapterClip, clipMode=clipMode, ...)    
}

"sread" <- function(object, clipMode,...){
	if (inherits(object,"SffReads")){
		if (missing(clipMode)) { clipMode <- clipMode(object) }
		if(!(clipMode %in% c("Full","Quality","Raw"))) error("clipMode must be one of Full, Quality, Raw")
		clipFull <- function(object){
			clipL <- pmax(1, pmax(start(qualityClip(object)),start(adapterClip(object)) )) 
			clipR <- pmin( 
				ifelse(end(qualityClip(object)) == 0 , width(object@sread),end(qualityClip(object))), 
				ifelse(end(adapterClip(object)) == 0 , width(object@sread), end(adapterClip(object))) )
			clipR <- pmax(clipL,clipR) #
			subseq(object@sread,start=clipL,end=clipR)
		}
		clipQuality <- function(object){
			clipL <- pmax(1, pmax(start(qualityClip(object)) )) 
			clipR <- ifelse(end(qualityClip(object)) == 0 , width(object),end(qualityClip(object)))
			clipR <- pmax(clipL,clipR)
			subseq(object@sread,start=clipL,end=clipR)
		}
		switch(clipMode,
		    "Full"=clipFull(object),
		    "Quality"=clipQuality(object),
		    "Raw"=object@sread)
	} else object@sread
}

setMethod(writeFasta, "SffReads",
          function(object, file, ...)
{
    dna <- sread(object)
    callGeneric(dna, file=file, ...)
})

### Accessor functions

setMethod(id, "SffReads", function(object) BStringSet(names(object@sread)))

setMethod(length, "SffReads", function(x) length(sread(x)))

setMethod(width, "SffReads", function(x) width(sread(x)))

setMethod(adapterClip, "SffReads", function(object) object@adapterClip)

setReplaceMethod( f="adapterClip",signature="SffReads", 
    definition=function(object,value){
	    if (class(value) != "IRanges")
			error("value must be of type IRanges object")
        object@adapterClip <-value 
        return (object)
})

setMethod(qualityClip, "SffReads", function(object) object@qualityClip)

setReplaceMethod( f="qualityClip",signature="SffReads", 
    definition=function(object,value){
	    if (class(value) != "IRanges")
			error("value must be of type IRanges object")
        object@qualityClip <-value 
        return (object)
})

setMethod(clipMode, "SffReads", function(object) object@clipMode)

setReplaceMethod( f="clipMode",signature="SffReads", 
    definition=function(object,value){
	    if (!(value %in% c("Full","Quality","Raw")))
			error("Unknown clipMode, must be one of 'Full','Quality','Raw'")
        object@clipMode <-value 
        return (object)
})


## coerce

setMethod(pairwiseAlignment, "SffReads",
          function(pattern, subject, ...)
          {
            pairwiseAlignment(sread(pattern), subject, ...)
})

## subset

setMethod("[", c("SffReads", "missing", "missing"),
          function(x, i, j, ..., drop=NA) 
			error("UserSubset:'[' must be called with only subscript 'i'")
)

setMethod("[", c("SffReads", "missing", "ANY"),
          function(x, i, j, ..., drop=NA)
			error("UserSubset:'[' must be called with only subscript 'i'")
)

setMethod("[", c("SffReads", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA)
			error("UserSubset:'[' must be called with only subscript 'i'")
)

.SffReads_subset <- function(x, i, j, ..., drop=TRUE) {
    if (length(list(...)) != 0L) 
		error("UserSubset:'[' must be called with only subscript 'i'")
	initialize(x, sread=sread(x)[i],
	               qualityClip=qualityClip(x)[i],
	               adapterClip=adapterClip(x)[i],
	##TODO:subset header
	               header=header(x),clipMode=clipMode(x))
}

setMethod("[", c(x="SffReads", i="ANY", j="missing"),
          .SffReads_subset)

setMethod(append, c("SffReads", "SffReads", "missing"),
    function(x, values, after=length(x)) 
{
	    initialize(x,
	               sread=append(sread(x), sread(values)),
				   qualityClip=append(qualityClip(x),qualityClip(values)),
				   adapterClip=append(adapterClip(x),adapterClip(values)),
	##TODO:add append headers to methods_SffHeader
	               header=list(header(x),header(value)),clipMode=clipMode(x))
})

setMethod(alphabetByCycle, "SffReads", ShortRead:::.abc_ShortRead)

setMethod(clean, "SffReads", function(object, ...) {
    alf <- alphabetFrequency(sread(object), baseOnly=TRUE)
    object[alf[,'other'] == 0]
})

setMethod(dustyScore, "SffReads", function(x, batchSize=NA, ...) {
    callGeneric(sread(x), batchSize=batchSize, ...)
})

setMethod(srorder, "SffReads", function(x, ...) callGeneric(sread(x), ...))

setMethod(srrank, "SffReads",  function(x, ...) callGeneric(sread(x), ...))

setMethod(srduplicated, "SffReads", function(x, ...) callGeneric(sread(x), ...))

setMethod(srdistance, c("SffReads", "ANY"), function(pattern, subject, ...){
		callGeneric(sread(pattern), subject, ...)
})

setMethod(srsort, "SffReads", function(x, ...) {
		x[srorder(x, ...)]
})


setMethod(tables, "SffReads", function(x, n=50, ...) {
		callGeneric(sread(x), n=n, ...)
})


###TODO: Fix qualityClip and adapterClip Narrow
setMethod(narrow, "SffReads",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
{
    initialize(x,
               sread=narrow(sread(x), start, end, width, use.names),
               qualityClip=narrow(qualityClip(x),start,end,width,use.names),
               adapterClip=narrow(adapterClip(x),start,end,width,use.names),
               header=header(x),clipMode=clipMode(x))
})


setMethod(trimLRPatterns, c(subject="SffReads"),
    function (Lpattern = "", Rpattern = "", subject, max.Lmismatch =
              0, max.Rmismatch = 0, with.Lindels = FALSE, with.Rindels
              = FALSE, Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
{
    ret <-
        callGeneric(Lpattern, Rpattern, sread(subject), max.Lmismatch,
                    max.Rmismatch, with.Lindels, with.Rindels, Lfixed,
                    Rfixed, ranges=TRUE)
    if (ranges)
        ret
    else 
        narrow(subject, start(ret), end(ret))
})

setMethod(trimEnds, "SffReads",
    function(object, a, left=TRUE, right=TRUE, relation=c("<=", "=="),
             ..., ranges=FALSE)
{
    rng <- callGeneric(sread(object), a, left, right, relation,
                       ..., ranges=TRUE)
    if (ranges) rng
    else narrow(object, start(rng), end(rng))
})

## manip
## show
setMethod(show, "SffReads", function(object) {
    callNextMethod()
    wd <- sort(unique(width(object)))
    if (length(wd)>2) wd <- paste(range(wd), collapse="..")
    cat("length:", length(object), "reads; width:", wd, "basepair; clipping mode:", clipMode(object), "\n")
})

## detail
setMethod(detail, "SffReads", function(x, ...) {
    callNextMethod()
    cat("\nClip Mode: ",clipMode(x),"\n")
    cat("\nsread:\n")
    show(sread(x))
    cat("\nclip points:\n")
    show(qualityClip(x))
    show(adapterClip(x))
})
