
## Inspector
setMethod(.sffValidity,
          signature=c("SffReads"),
          function(object) {
            #	cat("## SFFreads Validity ###\n")
            msg <- NULL
            object@qualityIR <- .solveIRangeSEW(width(object@sread),object@qualityIR)
            object@adapterIR <- .solveIRangeSEW(width(object@sread), object@adapterIR)
            object@customIR <- .solveIRangeSEW(width(object@sread),object@customIR)
            if (!(object@clipMode %in% availableClipModes(object))){
              txt <- sprintf(paste("wrong mode type must be one of",paste(availableClipModes(object),collapse=" ")))
              msg <- c(msg,txt)
            }
            cat("called SffReads validity")
            if (is.null(msg)) TRUE else msg
})

## constructor, missing customIR for now
"SffReads" <- 
function(sread, qualityIR, adapterIR, customIR, clipMode="raw", header)
{
  if (missing(header)) header = list()
  if (missing(qualityIR)) qualityIR=IRanges()
  if (missing(adapterIR)) adapterIR=IRanges()
  if (missing(customIR)) customIR=IRanges()
  
  new("SffReads", header=header, sread=sread,
      qualityIR=qualityIR, adapterIR= adapterIR, customIR=customIR, clipMode=clipMode)    
}

### Accessor functions

"sread" <- 
  function(object, start=NULL,end=NULL,width=NULL,clipmode,...){
	if (inherits(object,"SffReads")){
    IR <- solveSffSEW(object,start,end,width,clipmode,...)
    subseq(object@sread,start=start(IR),end=end(IR))
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
setMethod(names, 
          signature="SffReads", 
          function(x) names(x@sread))

### Reassign Read Names
setReplaceMethod( f="names",
                  signature="SffReads",
                  function(x,value){names(x@sread) <- value; return(x)})

### Print out Read Names as BString Set (ShortRead way)
setMethod(id,
          signature="SffReads",
          function(object){ 
            if (is.null(names(object@sread))) return(NULL)
            BStringSet(names(object@sread))
})

### number of sequences
setMethod(length, "SffReads", function(x) length(x@sread))

### vector of widths, is dependant on current clipMode
setMethod(width, "SffReads", function(x) width(solveSffSEW(x)))

##### Clipping Points
### IRange object of adapterClip points
setMethod(adapterClip, signature="SffReads", definition=function(object) .solveIRangeSEW(width(object@sread),object@adapterIR))

### Reassign adapterClip points
setReplaceMethod( f="adapterClip",signature="SffReads", 
                  definition=function(object,value){
                    if (class(value) != "IRanges")
                      stop("value must be of type IRanges object")
                    object@adapterIR <- .solveIRangeSEW(width(object@sread),value)
                    return (object)
})

### IRange object of qualityClip points
setMethod(qualityClip, "SffReads", function(object) .solveIRangeSEW(width(object@sread),object@qualityIR))
  
### Reassign qualityClip Points
setReplaceMethod( f="qualityClip",signature="SffReads", 
                  definition=function(object,value){
                    if (class(value) != "IRanges")
                      stop("value must be of type IRanges object")
                    object@qualityIR <- .solveIRangeSEW(width(object@sread),value) 
                    return (object)
                  })


### IRange object of adapterClip points
setMethod(customClip, "SffReads", function(object) .solveIRangeSEW(width(object@sread),object@customIR))

### Assign customClip points
setReplaceMethod( f="customClip",signature="SffReads", 
                  definition=function(object,value){
                    if (class(value) != "IRanges")
                      stop("value must be of type IRanges object")
                    ##TODO: Need to check validity of new clip points
                    object@customIR <- .solveIRangeSEW(width(object@sread),value)
                    return (object)
                  })

setMethod(rawClip,"SffReads", function(object) solveUserSEW(width(object@sread)))
setMethod(fullClip,"SffReads", function(object) solveSffSEW(object,clipMode="full"))

### Get the current clipMode for the object
setMethod(clipMode, "SffReads", function(object) object@clipMode)

### Reset the clip mode
setReplaceMethod( f="clipMode",signature="SffReads", 
                  definition=function(object,value){
                    if (!(value %in% availableClipModes(object)))
                      stop("clipMode not available, must be one of ",paste(availableClipModes(object),collapse=","))
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
	               qualityIR=qualityClip(x)[i],
	               adapterIR=adapterClip(x)[i],
	               customIR=if(length(x@customIR) != 0){customClip(x)[i]}else{IRanges()},
	             ##TODO:subset header
	               header=header(x),clipMode=clipMode(x))
}

setMethod("[", c(x="SffReads", i="ANY", j="missing"),
          .SffReads_subset)

setMethod(append, c("SffReads", "SffReads", "missing"),
    function(x, values, after=length(x)) 
{
      appendCustom <- function(IR1,IR2,width1,width2){
        if (length(IR1) != 0 & length(IR2) != 0)
          append(IR1,IR2)
        else if (length(IR1) != 0 & length(IR2)==0)
          append(IR1,IRanges(1,width2))
        else if (length(IR1) == 0 & length(IR2 !=0))
          append(IRanges(1,width1),IR2)
        else IRanges()      
      }
      
	    initialize(x,
           sread=append(x@sread, values@sread),
				   qualityIR=append(qualityClip(x),qualityClip(values)),
				   adapterIR=append(adapterClip(x),adapterClip(values)),
	         customIR=appendCustom(x@customIR,values@customIR,width(x),width(values)),
	         ##TODO:add append headers to methods_SffHeader
	         header=list(header(x),header(values)),clipMode=clipMode(x))
})

#### FIX
# 
# setMethod(reverseComplement, "SffReads",function(x, index, ...)
# {
#   if (missing(index)) index <- seq.int(1L,length(object))
#   if (is.logical(index)) index <- which(index)
#   if (!is.numeric(index)) stop("index must be either missing, a logical vector, or numeric vector")
#   newsff <- x
#   newsff@sread[index] <- reverseComplement(newsff@sread[index])
#   qualityClip(newsff)[index] <- IRanges(end=width(newsff@sread[index]) - start(qualityClip(newsff)[index])+1,
#                                         start  =width(newsff@sread[index]) - end(qualityClip(newsff)[index])+1)
#   adapterClip(newsff)[index] <- IRanges(end=width(newsff@sread[index]) - start(adapterClip(newsff)[index])+1,
#                                         start  =width(newsff@sread[index]) - end(adapterClip(newsff)[index])+1)
#   newsff
# })

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

setMethod(writeFasta, c("SffReads"),
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
