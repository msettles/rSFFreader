
### Accessor functions
setMethod(header, 
	        signature="SffHeader", 
	        function(object) object@header
)

## show
setMethod(show, 
          signature="SffHeader",
	        function(object) {
	          cat("class: ", class(object), "\n", sep="")
#           sapply(header(object), function(x) cat("file:",x$filename,"; number of reads:",x$number_of_reads,"; numer of flows:",x$number_of_flows,"\n"))
})

## detail
setMethod(detail,
	        signature="SffHeader",
	        function(x, ...) {
	          cat("class: ", class(x), "\n", sep="")
            cat("\nheader information:\n")
            show(header(x))   
})

