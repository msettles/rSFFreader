
setGeneric("header", function(object, ...)
           standardGeneric("header"))

setGeneric("adapterClip", function(object, ...)
           standardGeneric("adapterClip"))

setGeneric("qualityClip", function(object, ...)
           standardGeneric("qualityClip"))

setGeneric("clipMode", function(object, ...)
           standardGeneric("clipMode"))

setGeneric("writePhredQual", function(object, file, mode="w", ...)
           standardGeneric("writePhredQual"))

setGeneric("writeFastaQual", function(object, file, mode="w", ...)
           standardGeneric("writeFastaQual"))

setGeneric("clipMode<-",function(object,value) {standardGeneric("clipMode<-")} )
