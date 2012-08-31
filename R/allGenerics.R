
setGeneric("header", function(object, ...)
           standardGeneric("header"))

#### clipModes 
setGeneric("fullClip",function(object, ...)
  standardGeneric("fullClip"))
### does not get a replace method

setGeneric("adapterClip", function(object, ...)
  standardGeneric("adapterClip"))
setGeneric("adapterClip<-",function(object,value) {standardGeneric("adapterClip<-")} )

setGeneric("qualityClip", function(object, ...)
  standardGeneric("qualityClip"))
setGeneric("qualityClip<-",function(object,value) {standardGeneric("qualityClip<-")} )

setGeneric("rawClip", function(object, ...)
  standardGeneric("rawClip"))
### does not get a replace method

setGeneric("customClip", function(object, ...)
  standardGeneric("customClip"))
setGeneric("customClip<-",function(object,value) {standardGeneric("customClip<-")} )

setGeneric("clipMode", function(object, ...)
  standardGeneric("clipMode"))
setGeneric("clipMode<-",function(object,value) {standardGeneric("clipMode<-")} )

#setGeneric("id<-",function(x,value) {standardGeneric("id<-")} )
#setGeneric("sread<-",function(object,value) {standardGeneric("sread<-")} )
#setGeneric("quality<-",function(object,value) {standardGeneric("quality<-")} )

setGeneric("writePhredQual", function(object, filepath, mode="w", ...)
  standardGeneric("writePhredQual"))

## breaks when a file already exists, doesn't overwrite
setGeneric("writeFastaQual", function(object, basefilename, append=FALSE, ...)
  standardGeneric("writeFastaQual"))
