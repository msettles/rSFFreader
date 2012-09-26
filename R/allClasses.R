.sffValidity <- function(object) TRUE
setGeneric(".sffValidity")

## Classes
setClass("SffHeader", contains=c(),
         representation=representation(
           header="list"),
         prototype=prototype(
           header=list())
)

setClass("SffReads", contains="SffHeader",
         representation = representation(
	         sread="DNAStringSet",
           qualityIR="IRanges",
           adapterIR="IRanges",
           customIR="IRanges",
           clipMode="character"),
         prototype=prototype(
	       sread=DNAStringSet(character(0)),
           qualityIR=IRanges(numeric(0)),
           adapterIR=IRanges(numeric(0)),
           customIR=IRanges(numeric(0)),
           clipMode="Full"),
         validity=.sffValidity
)

setClass("SffReadsQ", contains="SffReads",
         representation=representation(
           quality="QualityScore"),
         prototype=prototype(
           quality=NumericQuality()),
         validity=.sffValidity
)
