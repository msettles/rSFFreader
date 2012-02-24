
### Accessor functions
setMethod(header, 
	      signature="SffHeader", 
	      function(object) object@header
)

## show
setMethod(show,
	      signature="SffHeader",
	      function(object) {
            callNextMethod()
            sapply(header(object), function(x) cat("file:",x$filename,"; number of reads:",x$number_of_reads,"; numer of flows:",x$number_of_flows,"\n"))
})

## detail
setMethod(detail,
	      signature="SffHeader",
	      function(x, ...) {
			callNextMethod()
			cat("\nheader information:\n")
    		show(header(x))   
})

