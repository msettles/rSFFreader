
#dyn.load("Roche454Reads/src/RocheSFF-io.so")
#.Call("read_roche_sff","Roche454Reads/inst/test/SmallTest.sff",as.integer(TRUE),"Roche454Reads")

readsff <- function(filenames, use.names=TRUE, verbose=TRUE){

    stopifnot(file.exists(filenames))
	if (!isTRUEorFALSE(use.names))
		stop("'use.names' must be TRUE or FALSE")
	if (!isTRUEorFALSE(verbose))
		stop("'verbose' must be TRUE or FALSE")
	
		lkup_seq <- get_xsbasetypes_conversion_lookup("B", "DNA")
	ans <- .Call("read_sff",filenames,use.names,lkup_seq, NULL,verbose,"rSFFreader")
	SffReads(ans[[2]],ans[[3]],ans[[4]],ans[[5]],"Full")
}

readsffheader <- function(filenames,verbose=TRUE) {
    stopifnot(file.exists(filenames))
	if (!isTRUEorFALSE(verbose))
		stop("'verbose' must be TRUE or FALSE")
	
	return(.Call("read_sff_header", filenames,verbose,"rSFFreader"))
}

readsffgeometry <- function(filenames) {
	stopifnot(file.exists(filenames))
	return(.Call("read_geometry", filenames,"rSFFreader"))
}
