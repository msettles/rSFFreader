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
    cat("called SffReadsQ validity")
    
    if (is.null(msg)) TRUE else msg
})


## constructor
SffReadsQ <- 
function(sread, quality, qualityIR, adapterIR, customIR, clipMode="raw", header)
{
    if (missing(header)) header = list()
    if (missing(qualityIR)) qualityIR=IRanges()
    if (missing(adapterIR)) adapterIR=IRanges()
    if (missing(customIR)) customIR=IRanges()

    if (class(quality) == "BStringSet") 
        quality <- FastqQuality(quality)
    if (class(quality) != "FastqQuality") stop("quality slot must be of type FastqQuality or BStringSet")

    new("SffReadsQ", header=header, sread=sread, quality=quality, 
        qualityIR=qualityIR, adapterIR= adapterIR, customIR=customIR, clipMode=clipMode)    
}

"quality" <- function(object, start=NULL, end=NULL, width=NULL, clipmode, ...){
  if (inherits(object,"SffReadsQ")){
    IR <- solveSffSEW(object,start,end,width,clipmode,...)
    FastqQuality(subseq(quality(object@quality),start=start(IR),end=end(IR)))
  } else object@quality
}
# 
# setReplaceMethod( f="quality",signature="SffReads", 
#     definition=function(object,value){
# 	    if (class(value) == "BStringSet") value <- FastqQuality(value)
# 	    if (class(value) != "FastqQuality")
# 			stop("value must be of type BStringSet or FastqQuality object")
# #TODO: More Checks on IRanges, length, names, etc.
#         object@quality <-value 
#         return (object)
# })
###### FIX
# setMethod(reverseComplement, "SffReadsQ",function(x, index, ...)
# {
#     if (missing(index)) index <- seq.int(1L,length(object))
# 		if (is.logical(index)) index <- which(index)
# 		if (!is.numeric(index)) stop("index must be either missing, a logical vector, or numeric vector")
# 		newsff <- x
# 		newsff@sread[index] <- reverseComplement(newsff@sread[index])
# 		qual <- quality(newsff@quality)
# 		qual[index] <- reverse(quality(newsff@quality)[index])
# 		newsff@quality <- FastqQuality(qual)
# 		qualityClip(newsff)[index] <- IRanges(end=width(newsff@sread[index]) - start(qualityClip(newsff)[index])+1,
# 												 start  =width(newsff@sread[index]) - end(qualityClip(newsff)[index])+1)
# 		adapterClip(newsff)[index] <- IRanges(end=width(newsff@sread[index]) - start(adapterClip(newsff)[index])+1,
# 												 start  =width(newsff@sread[index]) - end(adapterClip(newsff)[index])+1)
# 		newsff
# })
# 

#setMethod(pairwiseAlignment, "SffReadsQ",
#          function(pattern, subject, ...)
#          {
#            mc <- as.list(match.call())
#            if (is.null(mc$patternQuality))
#               mc$patternQuality <- quality(quality(pattern))
#             do.call(callNextMethod, c(list(pattern, subject), mc))
#           })

## subset

setMethod("[", c("SffReadsQ", "missing", "missing"),
          function(x, i, j, ..., drop=NA) 
	    stop("UserSubset:'[' must be called with only subscript 'i'")
)

setMethod("[", c("SffReadsQ", "missing", "ANY"),
          function(x, i, j, ..., drop=NA)
		stop("UserSubset:'[' must be called with only subscript 'i'")
)

setMethod("[", c("SffReadsQ", "ANY", "ANY"),
          function(x, i, j, ..., drop=NA)
		stop("UserSubset:'[' must be called with only subscript 'i'")
)


.SffReadsQ_subset <- function(x, i, j, ..., drop=TRUE) {
    if (0L != length(list(...)))
		stop("UserSubset:'[' must be called with only subscript 'i'")
    initialize(x, sread=x@sread[i],
               quality=x@quality[i],
               qualityIR=qualityClip(x)[i],
               adapterIR=adapterClip(x)[i],
               customIR=if(length(x@customIR) != 0){customClip(x)[i]}else{IRanges()},
##TODO:subset header
               header=header(x),clipMode=clipMode(x))
}

setMethod("[", c("SffReadsQ", "ANY", "missing"), .SffReadsQ_subset)

setMethod(append, c("SffReadsQ", "SffReadsQ", "missing"),
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
               quality=append(x@quality, values@quality),
			         qualityIR=append(qualityClip(x),qualityClip(values)),
			         adapterIR=append(adapterClip(x),adapterClip(values)),
               customIR=appendCustom(x@customIR,values@customIR,width(x),width(values)),
##TODO:add append headers to methods_SffHeader
               header=list(header(x),header(values)),clipMode=clipMode(x))
			
})

#setMethod(alphabetByCycle, "SffReadsQ", ShortRead::.abc_ShortReadQ)

#setMethod(alphabetScore, "SffReadsQ", ShortRead:::.forward_objq)


###TODO: Fix qualityClip and adapterClip Narrow
#setMethod(narrow, "SffReadsQ",
#    function(x, start=NA, end=NA, width=NA, use.names=TRUE)
#{
#    initialize(x,
#               sread=narrow(sread(x), start, end, width, use.names),
#               quality=narrow(quality(x), start, end, width, use.names),
#               qualityClip=qualityClip(x),
#               adapterClip=adapterClip(x),
#               header=header(x),clipMode=clipMode(x))
#})

#setMethod(trimTailw, "SffReadsQ",
#    function(object, k, a, halfwidth, ..., ranges=FALSE)
#{

#    rng <- callGeneric(quality(object), k, a, halfwidth, ...,
#                       ranges=TRUE)
#    if (ranges) rng
#    else narrow(object, 1L, end(rng))[0L != width(rng)]
#})

#setMethod(trimTails, "SffReadsQ",
#    function(object, k, a, successive=FALSE, ..., ranges=FALSE)
#{
#    rng <- callGeneric(quality(object), k, a, successive,
#                       ..., ranges=TRUE)
#    if (ranges) rng
#    else narrow(object, 1L, end(rng))[0L != width(rng)]
#})

#setMethod(trimEnds, "SffReadsQ",
#    function(object, a, left=TRUE, right=TRUE, relation=c("<=", "=="),
#             ..., ranges=FALSE)
#{
#    rng <- callGeneric(quality(object), a, left, right, relation,
#                       ..., ranges=TRUE)
#    if (ranges) rng
#    else narrow(object, start(rng), end(rng))
#})

### Functions to write out data
setMethod(writeFastq, "SffReadsQ", function(object, file, mode="w", full=FALSE, ...) {
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
          quality(quality(object)), file, mode,full, max_width)
    invisible(length(object))
})

setMethod(writePhredQual, "SffReadsQ", function(object, filepath, mode="w", ...) {
    if (length(filepath) != 1)
        sprintf("UserArgumentMismatch:'%s' must be '%s'",
                       "file", "character(1)")
    file <- path.expand(filepath)
    if (file.exists(file) && mode != "a")
        sprintf("UserArgumentMismatch:file '%s' exists, but mode is not 'a'",filepath)
    ## FIXME: different quality types
    max_width <- max( unique(width(quality(object))))
    
    .Call("write_phred_quality", id(object), 
          quality(quality(object)), file, mode, max_width)
    invisible(length(object))
})

setMethod(writeFastaQual, "SffReadsQ", function(object, basefilename, append=FALSE,  ...) {
    if (length(basefilename) != 1)
        sprintf("UserArgumentMismatch:'%s' must be '%s'",
                       "file", "character(1)")
    if (append == FALSE){mode="w"}else{mode="a"}
    file <- path.expand(basefilename)
    if ((file.exists(paste(basefilename,"fasta",sep=".")) && mode != "a") |
        (file.exists(paste(basefilename,"fasta.qual",sep=".")) && mode != "a"))
        sprintf("Warning:file '%s' exists, but mode is not 'a'", file)
    ## FIXME: different quality types ??
    max_width <- max(c(unique(width(names(sread(object)))),
                       unique(width(sread(object))),
                       unique(width(quality(object)))))
    writeXStringSet(sread(object), paste(file,"fasta",sep="."),append,format="fasta",...) ##
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
