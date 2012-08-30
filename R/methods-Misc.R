availableClipModes <- function(x) c("Full","Quality","Raw","Custom")

.solveIRangeSEW <- function(widths,IR){
  solveUserSEW(widths,start(IR),end(IR))
}

## Determine view on the data, called anytime one views the reads or qualitys
"solveSffSEW" <- function(object,starts=NULL,ends=NULL,widths=NULL,clipMode)
{
  if (!is.null(starts) || !is.null(ends) || !is.null(widths)) { 
    clipmode <- "Specified"
  } else if (!missing(clipMode)) {
    if (!(clipMode %in% availableClipModes()))
      txt <- sprintf("wrong mode type must be one of Full, Quality, Raw, or Custom")
    clipmode <- clipMode
  } else clipmode <- object@clipMode
  
  clipFull <- function(object){
    clipL <- pmax(start(qualityClip(object)),start(adapterClip(object)))
    clipR <- pmin(end(qualityClip(object)),end(adapterClip(object)))
    return(solveUserSEW(width(object@sread),clipL,clipR))
  }
  switch(clipmode,
         "Specified" = solveUserSEW(width(object@sread),starts,ends,widths),
         "Custom" = .solveIRangeSEW(width(object@sread),customClip(object)),
         "Full"=clipFull(object),
         "Quality"=.solveIRangeSEW(width(object@sread),qualityClip(object)),
         "Raw"=solveUserSEW(width(object@sread)))
}


setMethod(writePhredQual, "FastqQuality", function(object, file, mode="w", ...) {
  if (length(file) != 1)
    sprintf("UserArgumentMismatch:'%s' must be '%s'",
            "file", "character(1)")
  if (file.exists(file) && mode != "a")
    sprintf("UserArgumentMismatch:file '%s' exists, but mode is not 'a'",
            file)
  file <- path.expand(file)
  ## FIXME: different quality types
  max_width <- max(c(unique(width(names(sread(object)))),
                     unique(width(quality(object)))))
  .Call("write_phred_quality", id(object), 
        quality(object), file, mode, max_width)
  invisible(length(object))
})
