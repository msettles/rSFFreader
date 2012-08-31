
## Inspector
setMethod(.sffValidity, "SffReads", function(object) {
#	cat("## SFFreads Validity ###\n")
    msg <- NULL
    lens <- length(object@sread)
    lenqc <- length(object@qualityIR)
    lenac <- length(object@adapterIR)
    if (length(unique(c(lens,lenqc,lenac))) != 1) {
        txt <- sprintf("mismatch length in sread, qualityClip, or adapterClip: %d %d %d",
                       lens,lenqc,lenac)
        msg <- c(msg, txt)
    }
    
##TODO: need to add check for IRanges out of bounds
	if (!(object@clipMode %in% availableClipModes())){
		txt <- sprintf("wrong mode type must be one of",paste(availableClipModes(),collapse=" "))
		msg <- c(msg,txt)
	}
    if (is.null(msg)) TRUE else msg
})

## constructor
"SffReads" <- function(sread, qualityClip, adapterClip,
	clipMode=availableClipModes(), header, ...)
{
    if (missing(header)) header = list()
    clipMode = match.arg(clipMode)
    if (missing(qualityClip) | missing(adapterClip))
        emptyIR <- IRanges(start=rep(1,length(sread)),width=width(sread))
    if (missing(qualityClip)) qualityClip=emptyIR
    if (missing(adapterClip)) adapterClipj=emtpyIR
    
    new("SffReads", header=header, sread=sread,
        qualityClip=.solveIRangeSEW(width(sread),qualityClip),adapterClip= .solveIRangeSEW(width(sread), adapterClip), clipMode=clipMode, ...)    
}

### Accessor functions

"sread" <- function(object, clipmode,...){
	if (inherits(object,"SffReads")){
		if (missing(clipmode)) { clipmode <- clipMode(object) }
		if(!(clipmode %in% availableClipModes())) stop("clipmode must be one of",paste(availableClipModes(),collapse=" "))
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
		switch(clipmode,
		    "Full"=clipFull(object),
		    "Quality"=clipQuality(object),
		    "Raw"=object@sread)
	} else object@sread
}

# ### Replace sread object, must worry about clip points
# setReplaceMethod( f="sread",signature="SffReads", 
#                   definition=function(object,value){
#                     if (class(value) != "DNAStringSet")
#                       stop("value must be of type DNAStringSet object")
#                     #TODO: More Checks, Range objects
#                     object@sread <-value 
#                     return (object)
#                   })

### Print out Read Names
setMethod(names, "SffReads", function(x) names(x@sread))

### Reassign Read Names
setReplaceMethod( f="names",signature="SffReads",
    definition=function(x,value){names(x@sread) <- value; return(x)})

### Print out Read Names as BString Set (ShortRead way)
setMethod(id, "SffReads", function(object) BStringSet(names(object@sread)))


### Reassign Read Names
setReplaceMethod( f="id",signature="SffReads",
   definition=function(x,value){names(x@sread) <- value; return(x)})

### number of sequences
setMethod(length, "SffReads", function(x) length(x@sread))

### vector of widths, is dependant on current clipMode
setMethod(width, "SffReads", function(x) width(solveSffSEW(x)))

### IRange object of adapterClip points
setMethod(adapterClip, "SffReads", function(object) object@adapterIR)

### Reassign adapterClip points
setReplaceMethod( f="adapterClip",signature="SffReads", 
                  definition=function(object,value){
                    if (class(value) != "IRanges")
                      stop("value must be of type IRanges object")
                    object@adapterIR <- .solveIRangeSEW(width(object@sread),value)
                    return (object)
})

### IRange object of qualityClip points
setMethod(qualityClip, "SffReads", function(object) object@qualityIR)
  
### Reassign qualityClip Points
setReplaceMethod( f="qualityClip",signature="SffReads", 
                  definition=function(object,value){
                    if (class(value) != "IRanges")
                      stop("value must be of type IRanges object")
                    object@qualityIR <- .solveIRangeSEW(width(object@sread),value) 
                    return (object)
                  })


### IRange object of adapterClip points
setMethod(customClip, "SffReads", function(object) object@customIR)

### Reassign adapterClip points
setReplaceMethod( f="customClip",signature="SffReads", 
                  definition=function(object,value){
                    if (class(value) != "IRanges")
                      stop("value must be of type IRanges object")
                    ##TODO: Need to check validity of new clip points
                    object@customIR <- .solveIRangeSEW(width(object@sread),value)
                    return (object)
                  })

### Get the current clipMode for the object
setMethod(clipMode, "SffReads", function(object) object@clipMode)

### Reset the clip mode
setReplaceMethod( f="clipMode",signature="SffReads", 
                  definition=function(object,value){
                    if (!(value %in% availableClipModes()))
                      stop("clipMode must be one of ",paste(availableClipModes(),collapse=" "))
                    if (length(do.call(value,list(object))) == 0)
                      stop(paste("clipmMode:",value,"is not been set"))
                    object@clipMode <- value
                    return (object)
})

## subset

setMethod("[", c("SffReads", "missing", "missing"),
          function(x, i, j, ..., drop=NA) 
			stop("UserSubset:'[' must be called with only subscript 'i'")
)

setMethod("[", c("SffReads", "missing", "ANY"),
          function(x, i, j, ..., drop=NA)
			stop("UserSubset:'[' must be called with only subscript 'i'")
)

setMethod("[", c("SffReads", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA)
			stop("UserSubset:'[' must be called with only subscript 'i'")
)

.SffReads_subset <- function(x, i, j, ..., drop=TRUE) {
    if (length(list(...)) != 0L) 
		stop("UserSubset:'[' must be called with only subscript 'i'")
	  initialize(x, sread=x@sread[i],
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
           sread=append(x@sread, values@sread),
				   qualityClip=append(qualityClip(x),qualityClip(values)),
				   adapterClip=append(adapterClip(x),adapterClip(values)),
	##TODO:add append headers to methods_SffHeader
	         header=list(header(x),header(value)),clipMode=clipMode(x))
})

setMethod(reverseComplement, "SffReads",function(x, index, ...)
{
  if (missing(index)) index <- seq.int(1L,length(object))
  if (is.logical(index)) index <- which(index)
  if (!is.numeric(index)) stop("index must be either missing, a logical vector, or numeric vector")
  newsff <- x
  newsff@sread[index] <- reverseComplement(newsff@sread[index])
  qualityClip(newsff)[index] <- IRanges(end=width(newsff@sread[index]) - start(qualityClip(newsff)[index])+1,
                                        start  =width(newsff@sread[index]) - end(qualityClip(newsff)[index])+1)
  adapterClip(newsff)[index] <- IRanges(end=width(newsff@sread[index]) - start(adapterClip(newsff)[index])+1,
                                        start  =width(newsff@sread[index]) - end(adapterClip(newsff)[index])+1)
  newsff
})

##### Functions currently available in ShortRead that need implementation here

# setMethod(alphabetByCycle, "SffReads", ShortRead:::.abc_ShortRead)
# 
# setMethod(clean, "SffReads", function(object, ...) {
#     alf <- alphabetFrequency(sread(object), baseOnly=TRUE)
#     object[alf[,'other'] == 0]
# })
# 
# setMethod(dustyScore, "SffReads", function(x, batchSize=NA, ...) {
#     callGeneric(sread(x), batchSize=batchSize, ...)
# })
# 
# setMethod(srorder, "SffReads", function(x, ...) callGeneric(sread(x), ...))
# 
# setMethod(srrank, "SffReads",  function(x, ...) callGeneric(sread(x), ...))
# 
# setMethod(srduplicated, "SffReads", function(x, ...) callGeneric(sread(x), ...))
# 
# setMethod(srdistance, c("SffReads", "ANY"), function(pattern, subject, ...){
# 		callGeneric(sread(pattern), subject, ...)
# })
# 
# setMethod(srsort, "SffReads", function(x, ...) {
# 		x[srorder(x, ...)]
# })
# 
# setMethod(tables, "SffReads", function(x, n=50, ...) {
# 		callGeneric(sread(x), n=n, ...)
# })

# ## coerce
# setMethod(pairwiseAlignment, "SffReads",
#           function(pattern, subject, ...)
#           {
#             pairwiseAlignment(sread(pattern), subject, ...)
# })

###TODO: Fix qualityClip and adapterClip Narrow
#setMethod(narrow, "SffReads",
#    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
#{
#    initialize(x,
#               sread=narrow(sread(x), start, end, width, use.names),
#               qualityClip=narrow(qualityClip(x),start,end,width,use.names),
#               adapterClip=narrow(adapterClip(x),start,end,width,use.names),
#               header=header(x),clipMode=clipMode(x))
#})


#setMethod(trimLRPatterns, c(subject="SffReads"),
#    function (Lpattern = "", Rpattern = "", subject, max.Lmismatch =
#              0, max.Rmismatch = 0, with.Lindels = FALSE, with.Rindels
#              = FALSE, Lfixed = TRUE, Rfixed = TRUE, ranges = FALSE)
#{
#    ret <-
#        callGeneric(Lpattern, Rpattern, sread(subject), max.Lmismatch,
#                    max.Rmismatch, with.Lindels, with.Rindels, Lfixed,
#                    Rfixed, ranges=TRUE)
#    if (ranges)
#        ret
#    else 
#        narrow(subject, start(ret), end(ret))
#})

#setMethod(trimEnds, "SffReads",
#    function(object, a, left=TRUE, right=TRUE, relation=c("<=", "=="),
#             ..., ranges=FALSE)
#{
#    rng <- callGeneric(sread(object), a, left, right, relation,
#                       ..., ranges=TRUE)
#    if (ranges) rng
#    else narrow(object, start(rng), end(rng))
#})

setMethod(writeFasta, "SffReads",
          function(object, file, ...)
          {
            dna <- sread(object)
            callGeneric(dna, file=file, ...)
          })

## manip
## show
setMethod(show, "SffReads", function(object) {
    callNextMethod()
    wd <- paste(range(width(object)), collapse="..")
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
