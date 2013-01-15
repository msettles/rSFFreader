##Inspector for flowgrams
setMethod(.sffValidity, "SffReadsF", function(object) {
  msg <- NULL
  lenfg <- nrow(object@flowgram) #number of rows in the flowgram
  widfg <- ncol(object@flowgram) #number of columns in the flowgram
  lenind <- length(object@flowindices) #number of flow indices for current read
  widind <- sapply(object@flowindices, length) #vector of lenghts of each flow index per read
  
  #Four cases to check:
  # 1. Is the number of bases (in the header) equal to the number of indices?
  # 2. Is the number of flows (in the header) equal to the number of columns of the flowgram?
  # 3. Is the length of the seq object itself equal to the number of indices?
  # 4. Is the length of the seq object itself equal to the number of rows in the flowgram?
    
  #number of flow indices should be equal to the number of bases called in the reads
  if (all(width(object@sread) != widind)) {
    txt <- paste("Mimatch betwen the width of the flow indices and the width of the reads")
    msg <- c(msg, txt)
  }
  
  
  #number of columns should be equal to the number of flows
  if (header(object)[[1]]$number_of_flows_per_read != widfg) {
    txt <- paste("Mismatch between the number of flows in the header and number of columns in the flowgram:",
                 widfg, header(object)[[1]]$number_of_flows_per_read, sep=" ")
    msg <- c(msg, txt)
  }
  
  #The length of the sequence object should be equal to the length of the flow indices
  if (length(object@sread) != lenind) {
    txt <- paste("Mismatch between the length of the sequence data and the number of indices:",
               lenind, length(object@sread), sep=" ")
    msg <- c(msg, txt) 
  }
  
  if (length(object@sread) != lenfg) {
    txt <- paste("Mismatch between the length of the sequence data and the number of rows in the flowgram:",
                 lenfg, length(object@sread), sep=" ")
    msg <- c(msg, txt)
  }
  
#  #output
  if (is.null(msg)) TRUE
  else msg
}
  
          )


#constructor
SffReadsF <- 
  function(sread, quality, qualityIR, adapterIR, customIR, clipMode="raw", header, flowgram, flowindices)
  {
    if (missing(header)) header = list()
    if (missing(qualityIR)) qualityIR=IRanges()
    if (missing(adapterIR)) adapterIR=IRanges()
    if (missing(customIR)) customIR=IRanges()
    
    if (class(quality) == "BStringSet") 
      quality <- FastqQuality(quality)
    if (class(quality) != "FastqQuality") stop("quality slot must be of type FastqQuality or BStringSet")
    if (class(flowgram) != "matrix" && is.numeric(flowgram))
        stop("flowgram slot must be of type 'numeric matrix'")
    if (class(flowindices) != "list"&& is.numeric(unlist(flowindices)))
        stop("flowgram indices must be of type 'numeric list'")
    
    new("SffReadsF", header=header, sread=sread, quality=quality, 
        qualityIR=qualityIR, adapterIR= adapterIR, customIR=customIR, clipMode=clipMode, flowgram=flowgram, flowindices=flowindices)    
  }
## Print out flowgrams ###
setMethod(flowgram, 
          signature="SffReadsF", 
          function(object) object@flowgram)

### Print out flow indices ###
setMethod(flowindices, 
          signature="SffReadsF", 
          function(object) object@flowindices)

  