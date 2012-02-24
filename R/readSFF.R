
### debug
#lkup_seq <- get_xsbasetypes_conversion_lookup("B", "DNA")
#ans <- .Call("read_sff","rSFFreader/inst/test/SmallTest.sff",TRUE,lkup_seq, NULL,TRUE,"rSFFreader")
#new("SffReadsQ",sread=ans[[2]],quality=FastqQuality(ans[[3]]),
#  qualityClip=ans[[4]],adapterClip=ans[[5]],clipMode="Full",header=ans[[1]])

## functions 
## returns the contents of the SFF file into either a SffReads or SffReadsQ class, which acts and behaves similar to
## the ShortRead and ShortReadQ classes from package ShortRead
readsff <- function(filenames, use.qualities=TRUE, use.names=TRUE,clipMode=c("Full","Quality","Raw"), verbose=TRUE){

    if (!use.names) warning ("Currently use.names is not used, by default names will always be returned.")
    stopifnot(file.exists(filenames))
    clipMode <- match.arg(clipMode)
	if (!isTRUEorFALSE(use.names))
		stop("'use.names' must be TRUE or FALSE")
	if (!isTRUEorFALSE(verbose))
		stop("'verbose' must be TRUE or FALSE")
	
	lkup_seq <- get_xsbasetypes_conversion_lookup("B", "DNA")
	ans <- .Call("read_sff",filenames,use.names,lkup_seq, NULL,verbose,"rSFFreader")
	if (use.qualities){
    	SffReadsQ(sread=ans[["sread"]],quality=ans[["quality"]],
	    qualityClip=ans[["qualityClip"]],adapterClip=ans[["adapterClip"]],clipMode=clipMode,header=ans[["header"]])
    } else {
    	SffReads(ans[["sread"]],ans[["qualityClip"]],ans[["adapterClip"]],clipMode,ans[["header"]])
    }
}

## Returns a list of size 2
readsffgeometry <- function(filenames) {
	stopifnot(file.exists(filenames))
 	sffgeometry <- .Call("read_geometry", filenames,"rSFFreader")
    names(sffgeometry) <- c("nReads","Read_Widths")
	return(sffgeometry)
}


readsffheader <- function(filenames,verbose=TRUE) {
    stopifnot(file.exists(filenames))
    if (!isTRUEorFALSE(verbose))
		stop("'verbose' must be TRUE or FALSE")
	ans <- .Call("read_sff_header", filenames,verbose,"rSFFreader")
	new("SffHeader", header=ans)
}
