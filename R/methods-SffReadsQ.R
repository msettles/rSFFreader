
## Inspector
setMethod(.sffValidity, "SffReadsQ", function(object) {
    msg <- NULL
    lenq <- length(object@quality)
    lens <- length(object@sread)
    if (length(unique(c(lenq,lens))) != 1) {
        txt <- sprintf("mismatch length in quality and sread, quality: %d %d",
                       lens, lenq)
        msg <- c(msg, txt)
    }
    if (!all(width(object@quality) == width(object@sread))) {
        txt <- sprintf("some sread and quality widths differ")
        msg <- c(msg, txt)
    }
    if (is.null(msg)) TRUE else msg
})


## constructor
SffReadsQ <- function(sread, quality, qualityClip, adapterClip,
	clipMode=c("Full", "Quality", "Raw"), header, ...)
{
    if (missing(header)) header = list()
    clipMode = match.arg(clipMode)
    if (missing(qualityClip) | missing(adapterClip))
        emptyIR <- IRanges(start=integer(length(sread)),width=integer(length(sread)))
    if (missing(qualityClip)) qualityClip=emptyIR
    if (missing(adapterClip)) adapterClipj=emtpyIR
    
    if (class(quality) == "BStringSet") 
        quality <- FastqQuality(quality)
    if (class(quality) != "FastqQuality") error("quality slot must be of type FastqQuality or BStringSet")

    new("SffReadsQ", header=header, sread=sread, quality=quality, 
        qualityClip=qualityClip, adapterClip=adapterClip, clipMode=clipMode, ...)    
}

"quality" <- function(object, clipmode, ...){
	if (inherits(object,"SffReadsQ")){
		if (missing(clipmode)) { clipmode <- clipMode(object) }
		if(!(clipmode %in% c("Full","Quality","Raw"))) error("clipmode must be one of Full, Quality, Raw")
		clipFull <- function(object){
			clipL <- pmax(1, pmax(start(qualityClip(object)),start(adapterClip(object)) )) 
			clipR <- pmin( 
				ifelse(end(qualityClip(object)) == 0 , width(object@quality),end(qualityClip(object))), 
				ifelse(end(adapterClip(object)) == 0 , width(object@quality), end(adapterClip(object))) )
			FastqQuality(subseq(quality(object@quality),start=clipL,end=clipR))
		}
		clipQuality <- function(object){
			clipL <- pmax(1, pmax(start(qualityClip(object)) )) 
			clipR <- ifelse(end(qualityClip(object)) == 0 , width(object),end(qualityClip(object)))
			FastqQuality(subseq(quality(object@quality),start=clipL,end=clipR))
		}
		switch(clipmode,
		    "Full"=clipFull(object),
		    "Quality"=clipQuality(object),
		    "Raw"=object@quality)
	} else object@quality
}

setMethod(pairwiseAlignment, "SffReadsQ",
          function(pattern, subject, ...)
          {
            mc <- as.list(match.call())
            if (is.null(mc$patternQuality))
              mc$patternQuality <- quality(quality(pattern))
            do.call(callNextMethod, c(list(pattern, subject), mc))
          })

## subset

setMethod("[", c("SffReadsQ", "missing", "missing"),
          function(x, i, j, ..., drop=NA) 
	    error("UserSubset:'[' must be called with only subscript 'i'")
)

setMethod("[", c("SffReadsQ", "missing", "ANY"),
          function(x, i, j, ..., drop=NA)
		error("UserSubset:'[' must be called with only subscript 'i'")
)

setMethod("[", c("SffReadsQ", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA)
		error("UserSubset:'[' must be called with only subscript 'i'")
)

.SffReadsQ_subset <- function(x, i, j, ..., drop=TRUE) {
    if (0L != length(list(...)))
		error("UserSubset:'[' must be called with only subscript 'i'")
    initialize(x, sread=sread(x)[i],
               quality=quality(x)[i],
               qualityClip=qualityClip(x)[i],
               adapterClip=adapterClip(x)[i],
##TODO:subset header
               header=header(x),clipMode=clipMode(x))
}

setMethod("[", c("SffReadsQ", "ANY", "missing"), .SffReadsQ_subset)

setMethod(append, c("SffReadsQ", "SffReadsQ", "missing"),
    function(x, values, after=length(x))
{
    initialize(x,
               sread=append(sread(x), sread(values)),
               quality=append(quality(x), quality(values)),
			   qualityClip=append(qualityClip(x),qualityClip(values)),
			   adapterClip=append(adapterClip(x),adapterClip(values)),
##TODO:add append headers to methods_SffHeader
               header=list(header(x),header(value)),clipMode=clipMode(x))
			
})

###TODO: Fix qualityClip and adapterClip Narrow
setMethod(narrow, "SffReadsQ",
    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
{
    initialize(x,
               sread=narrow(sread(x), start, end, width, use.names),
               quality=narrow(quality(x), start, end, width, use.names),
               qualityClip=narrow(qualityClip(x),start,end,width,use.names),
               adapterClip=narrow(adapterClip(x),start,end,width,use.names),
               header=header(x),clipMode=clipMode(x))
})

setMethod(alphabetByCycle, "SffReadsQ", ShortRead:::.abc_ShortReadQ)

setMethod(alphabetScore, "SffReadsQ", ShortRead:::.forward_objq)

setMethod(trimTailw, "SffReadsQ",
    function(object, k, a, halfwidth, ..., ranges=FALSE)
{

    rng <- callGeneric(quality(object), k, a, halfwidth, ...,
                       ranges=TRUE)
    if (ranges) rng
    else narrow(object, 1L, end(rng))[0L != width(rng)]
})

setMethod(trimTails, "SffReadsQ",
    function(object, k, a, successive=FALSE, ..., ranges=FALSE)
{
    rng <- callGeneric(quality(object), k, a, successive,
                       ..., ranges=TRUE)
    if (ranges) rng
    else narrow(object, 1L, end(rng))[0L != width(rng)]
})

setMethod(trimEnds, "SffReadsQ",
    function(object, a, left=TRUE, right=TRUE, relation=c("<=", "=="),
             ..., ranges=FALSE)
{
    rng <- callGeneric(quality(object), a, left, right, relation,
                       ..., ranges=TRUE)
    if (ranges) rng
    else narrow(object, start(rng), end(rng))
})

### Functions to write out data

setMethod(writeFastq, "SffReadsQ", function(object, file, mode="w", ...) {
    if (length(file) != 1)
        sprintf("UserArgumentMismatch:'%s' must be '%s'",
                       "file", "character(1)")
    if (file.exists(file) && mode != "a")
        sprintf("UserArgumentMismatch:file '%s' exists, but mode is not 'a'",
                       file)
    file <- path.expand(file)
    ## FIXME: different quality types
    max_width <- max(c(unique(width(names(sread(object)))),
                       unique(width(sread(object))),
                       unique(width(quality(object)))))
    .Call(".write_fastq", id(object), sread(object),
          quality(quality(object)), file, mode, max_width)
    invisible(length(object))
})

setMethod(writePhredQual, "SffReadsQ", function(object, file, mode="w", ...) {
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
          quality(quality(object)), file, mode, max_width)
    invisible(length(object))
})

setMethod(writeFastaQual, "SffReadsQ", function(object, basefilename, mode="w", ...) {
    if (length(basefilename) != 1)
        sprintf("UserArgumentMismatch:'%s' must be '%s'",
                       "file", "character(1)")
    if ((file.exists(paste(basefilename,"fasta",sep=".")) && mode != "a") |
        (file.exists(paste(basefilename,"fasta.qual",sep=".")) && mode != "a"))
        sprintf("UserArgumentMismatch:file '%s' exists, but mode is not 'a'", file)
    file <- path.expand(basefilename)
    ## FIXME: different quality types
    max_width <- max(c(unique(width(names(sread(object)))),
                       unique(width(sread(object))),
                       unique(width(quality(object)))))
    write.XStringSet(sread(object), paste(file,"fasta",sep="."), ..., format="fasta")
    writePhredQual(object, paste(file,"fasta.qual",sep="."), mode, max_width)
    invisible(length(object))
})

### Write function "writeFastaQualXML" to produce sff_extract like output

## detail
setMethod(detail, "SffReadsQ", function(x, ...) {
	callNextMethod()
    cat("\nquality:\n")
    show(quality(x))    
})
