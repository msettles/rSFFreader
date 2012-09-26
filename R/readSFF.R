
### Sample Data
load454SampleData <- function() {
  readSff(system.file("extdata","Small454Test.sff",package="rSFFreader"))
}
loadIonSampleData <- function() {
  readSff(system.file("extdata","SmallTorrentTest.sff",package="rSFFreader"))
}

#lkup_seq <- get_seqtype_conversion_lookup("B", "DNA")
#ans <- .Call("read_sff","rSFFreader/inst/test/SmallTest.sff",TRUE,lkup_seq, NULL,TRUE,"rSFFreader")
#new("SffReadsQ",sread=ans[[2]],quality=FastqQuality(ans[[3]]),
#  qualityClip=ans[[4]],adapterClip=ans[[5]],clipMode="Full",header=ans[[1]])

### Some clip points (specifically adapter clip points) will have 0 values, reset to full length)
.fixSFFclipPoints2Iranges <- function(widths,IR){
  clipL <- pmax(1, start(IR)) 
  clipR <- ifelse(end(IR) == 0 , widths,end(IR))
  clipR <- pmin(widths,clipR)
  return(solveUserSEW(widths,clipL,clipR))  
}

## functions 
## returns the contents of the SFF file into either a SffReads or SffReadsQ class, which acts and behaves similar to
## the ShortRead and ShortReadQ classes from package ShortRead
readSff <- function(filenames, use.qualities=TRUE, use.names=TRUE,clipMode=c("full","adapter","quality","raw"), verbose=TRUE){
  if (!use.names) warning ("Currently use.names is not used, by default names will always be returned.")
  stopifnot(file.exists(filenames))
  clipMode <- match.arg(clipMode)
	if (!isTRUEorFALSE(use.names))
		stop("'use.names' must be TRUE or FALSE")
	if (!isTRUEorFALSE(verbose))
		stop("'verbose' must be TRUE or FALSE")
	lkup_seq <- get_seqtype_conversion_lookup("B", "DNA")
	ans <- .Call("read_sff",filenames,use.names,lkup_seq, NULL,verbose,"rSFFreader")
  widths <- width(ans[["sread"]])
	if (use.qualities){
    	SffReadsQ(sread=ans[["sread"]],quality=ans[["quality"]],
                qualityIR=.fixSFFclipPoints2Iranges(widths,ans[["qualityClip"]]),adapterIR=.fixSFFclipPoints2Iranges(widths,ans[["adapterClip"]]),
                clipMode=clipMode,header=ans[["header"]])
    } else {
    	SffReads(ans[["sread"]],
               qualityIR=.fixSFFclipPoints2Iranges(widths,ans[["qualityClip"]]),adapterIR=.fixSFFclipPoints2Iranges(widths,ans[["adapterClip"]]),
               clipMode=clipMode,ans[["header"]])
    }
}

## Returns a list of size 2
readSffGeometry <- function(filenames) {
	stopifnot(file.exists(filenames))
 	sffgeometry <- .Call("sff_geometry", filenames,"rSFFreader")
    names(sffgeometry) <- c("nReads","Read_Widths")
	return(sffgeometry)
}


readSffHeader <- function(filenames,verbose=TRUE) {
    stopifnot(file.exists(filenames))
    if (!isTRUEorFALSE(verbose))
		stop("'verbose' must be TRUE or FALSE")
	ans <- .Call("read_sff_header", filenames,verbose,"rSFFreader")
	new("SffHeader", header=ans)
}
