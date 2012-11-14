
setGeneric("header", function(object, ...)
           standardGeneric("header"))

#### clipModes 
setGeneric("rawClip", function(object, ...) standardGeneric("rawClip") )
### does not get a replace method

setGeneric("fullClip",function(object, ...) standardGeneric("fullClip") )
### does not get a replace method

setGeneric("adapterClip", function(object, ...) standardGeneric("adapterClip") )
setGeneric("adapterClip<-",function(object,value) standardGeneric("adapterClip<-") )

setGeneric("qualityClip", function(object, ...) standardGeneric("qualityClip") )
setGeneric("qualityClip<-",function(object,value) standardGeneric("qualityClip<-") )

setGeneric("customClip", function(object, ...) standardGeneric("customClip") )
setGeneric("customClip<-",function(object,value) standardGeneric("customClip<-") )

setGeneric("clipMode", function(object, ...) standardGeneric("clipMode") )
setGeneric("clipMode<-",function(object,value) standardGeneric("clipMode<-") )

setGeneric("writePhredQual", function(object, filepath, mode="w", ...)
  standardGeneric("writePhredQual"))

setGeneric("writeFastaQual", function(object, basefilename, append=FALSE, ...)
  standardGeneric("writeFastaQual"))
