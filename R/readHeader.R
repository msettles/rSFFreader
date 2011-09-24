
#dyn.load("Roche454Reads/src/RocheSFF-io.so")
#.Call("read_roche_sff","Roche454Reads/inst/test/SmallTest.sff",as.integer(TRUE),"Roche454Reads")

readsff <- function(filenames, use.names=TRUE, verbose=TRUE){
	lkup_seq <- as.integer(c( NA, NA, NA, NA, NA, NA, NA, NA, 
			       NA, NA, NA, NA, NA, NA, NA, NA, 
				   NA, NA, NA, NA, NA, NA, NA, NA, 
				   NA, NA, NA, NA, NA, NA, NA, NA, 
				   NA, NA, NA, NA, NA, NA, NA, NA, 
				   NA, NA, NA, 32, NA, 16, NA, NA,
				   NA, NA, NA, NA, NA, NA, NA, NA, 
				   NA, NA, NA, NA, NA, NA, NA, NA, 
				   NA, 1,  14, 2,  13, NA, NA, 4,  
				   11, NA, NA, 12, NA, 3,  15, NA, 
				   NA, NA, 5,  6,  8,  NA, 7,  9,  
				   NA, 10, NA, NA, NA, NA, NA, NA, 
				   NA, 1,  14, 2,  13, NA, NA, 4,  
				   11, NA, NA, 12, NA, 3,  15, NA, 
				   NA, NA, 5,  6,  8,  NA, 7,  9,  
				   NA, 10))
    stopifnot(file.exists(filenames))
	if (!isTRUEorFALSE(use.names))
		stop("'use.names' must be TRUE or FALSE")
	if (!isTRUEorFALSE(verbose))
		stop("'verbose' must be TRUE or FALSE")
	
	return(.Call("read_roche_sff",filenames,use.names,lkup_seq, NULL,verbose,"Roche454Reads"))
}

readSFFheader <- function(filenames,verbose=TRUE) {
    stopifnot(file.exists(filenames))
	if (!isTRUEorFALSE(verbose))
		stop("'verbose' must be TRUE or FALSE")
	
	return(.Call("read_roche_sff_header", filenames,verbose,"Roche454Reads"))
}
