

availableClipModes <- function() c("full","adapter","quality","raw","custom")

.solveIRangeSEW <- function(widths,IR){
  if (length(IR) == 0) stop("IRanges object of length 0")
  solveUserSEW(widths,start(IR),end(IR))
}

## Determine view on the data, called anytime one views the reads or qualitys
"solveSffSEW" <- function(object,starts=NULL,ends=NULL,widths=NULL,clipMode)
{
  if (!is.null(starts) || !is.null(ends) || !is.null(widths)) { 
    clipmode <- "Specified"
  } else if (!missing(clipMode)) {
    if (!(clipMode %in% availableClipModes()))
      txt <- sprintf(paste("wrong mode type must be one of",paste(availableClipModes(),collapse=",")))
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
  .Call("write_phred_quality", names(object), 
        quality(object), file, mode, max_width)
  invisible(length(object))
})
