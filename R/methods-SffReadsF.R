##Inspector for flowgrams
setMethod(.sffValidity, "SffReadsF", function(object) {
  msg <- NULL
  lenfg <- nrow(object@flowgram) #number of rows in the flowgram
  widfg <- ncol(object@flowgram) #number of columns in the flowgram
  lenind <- length(object@flowindicies) #number of flow indices for current read
  widind <- sapply(object@flowindicies, length) #vector of lenghts of each flow index per read
  
  #Four cases to check:
  # 1. Is the number of bases (in the header) equal to the number of indices?
  # 2. Is the number of flows (in the header) equal to the number of columns of the flowgram?
  # 3. Is the length of the seq object itself equal to the number of indices?
  # 4. Is the length of the seq object itself equal to the number of rows in the flowgram?
  
  
  
  
  #number of flow indices should be equal to the number of bases called
  if (header(object[[1]]$number_of_flows_per_read) != widind) {
    txt <- paste("Mimatch betwen the number of flows in the header and in the number of flow indices ")
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
               lenind, header(object)$number_of_flows_per_read, sep=" ")
    msg <- c(msg, txt) 
  }
  
#  ### add in above for lenfg
  #ind ! width(object@sread){
    
#  }
  
#  #output
  if (is.null(msg)) TRUE
  else msg
  
} 
          
)
  