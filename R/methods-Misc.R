

availableClipModes <- function(object){
  if (class(object) %in% c("SffReads","SffReadsQ")){
    avail <- c("raw")
    if (length(object@customIR)) avail <- c(avail,"custom")
    if (length(object@qualityIR)) avail <- c("adapter", avail)
    if (length(object@adapterIR)) avail <- c("quality", avail)
    if (length(object@qualityIR) & length(object@adapterIR)) avail <- c("full",avail)
    avail
    ## c("full","adapter","quality","raw","custom") 
  } else stop("clip modes only make sense for object of types SffReads or SffReadsQ classes")
}

.solveIRangeSEW <- function(widths,IR){
  if (length(IR) == 0) return(IRanges())
  solveUserSEW(widths,start(IR),end(IR))
}

## Determine view on the data, called anytime one views the reads or qualitys
"solveSffSEW" <- function(object,starts=NULL,ends=NULL,widths=NULL,clipMode)
{
  if (!is.null(starts) || !is.null(ends) || !is.null(widths)) { 
    clipmode <- "Specified"
  } else if (!missing(clipMode)) {
    if (!(clipMode %in% availableClipModes(object)))
      txt <- sprintf(paste("wrong mode type must be one of",paste(availableClipModes(object),collapse=",")))
    clipmode <- clipMode
  } else clipmode <- object@clipMode
  
  clipFull <- function(object){
    clipL <- pmax(start(object@qualityIR),start(object@adapterIR))
    clipR <- pmin(end(object@qualityIR),end(object@adapterIR))
    return(solveUserSEW(width(object@sread),clipL,clipR))
  }
  switch(clipmode,
         "specified" = solveUserSEW(width(object@sread),starts,ends,widths),
         "custom" = .solveIRangeSEW(width(object@sread),object@customIR),
         "full"=clipFull(object),
         "adapter"=.solveIRangeSEW(width(object@sread),object@adapterIR),
         "quality"=.solveIRangeSEW(width(object@sread),object@qualityIR),
         "raw"=solveUserSEW(width(object@sread)))
}

#### fix id, does not exist
setMethod(writePhredQual, "FastqQuality", function(object, filepath, mode="w", ...) {
  if (length(filepath) != 1)
    sprintf("UserArgumentMismatch:'%s' must be '%s'",
            "file", "character(1)")
  filepath <- path.expand(filepath)
  if (file.exists(filepath) && mode != "a")
    sprintf("UserArgumentMismatch:file '%s' exists, but mode is not 'a'",
            filepath)
  ## FIXME: different quality types
  max_width <- max( unique(width(quality(object))))
  .Call("write_phred_quality", names(object), quality(object), file, mode, max_width)
  invisible(length(object))
})
