
.sffValidity <- function(object) TRUE

setGeneric(".sffValidity")

setClass("SffHeader", 
         representation=representation(
           header="list"),
         prototype=prototype(
           header=list())
)


setClass("SffReads", #contains="SffHeader",
         representation=representation(
           sread="DNAStringSet",
           quality="QualityScore",
           qualityClip="IRanges",
           adapterClip="IRanges",
           clipMode="character"),
         prototype=prototype(
	       sread=DNAStringSet(character(0)),
           quality=NumericQuality(),
           qualityClip=IRanges(numeric(0)),
           adapterClip=IRanges(numeric(0)),
           clipMode="Full"),
         validity=.sffValidity
)
