
.sffValidity <- function(object) TRUE

setGeneric(".sffValidity")

## Virtual base classes

setClass(".SFFBase")

setClass("SffHeader", contains=".SFFBase",
         representation=representation(
           header="list"),
         prototype=prototype(
           header=list())
)

setClass("SffReads", contains="SffHeader",
         representation = representation(
	       sread="DNAStringSet",
           qualityClip="IRanges",
           adapterClip="IRanges",
           clipMode="character"),
         prototype=prototype(
	       sread=DNAStringSet(character(0)),
           qualityClip=IRanges(numeric(0)),
           adapterClip=IRanges(numeric(0)),
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

setClass("SffReadsF", contains="SffReadsQ",
         representation=representation(
           flowgram="matrix",
           flowindicies="list"),
         prototype=prototype(
           flowgram=matrix(numeric(0)),
           flowindicies=list(numeric(0))),
         validity=.sffValidity         
)
