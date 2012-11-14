\name{readSff}
\alias{readSff}

\title{Read SFF-formatted files into SffReads or SffReadsQ objects}

\description{
  \code{readSff} reads all sff files specified in \code{filenames}, returning an R object of type \code{\linkS4class{SffReads}} or 
  \code{\linkS4class{SffReadsQ}} which acts and behaves similarly to the \code{\linkS4class{ShortRead}} and \code{\linkS4class{ShortReadQ}} 
  classes from package ShortRead
}

\usage{
readSff(filenames, use.qualities=TRUE, use.names=TRUE, clipMode = c("full", "adapter", "quality", "raw"), verbose=TRUE)
}


\arguments{
  \item{filenames}{File or files to be read in.}
  
  \item{use.qualities}{\code{logical(1)} indicating whether to include quality values in the returned object.}
  
  \item{use.names}{\code{logical(1)} indicating whether to return read names, not yet implemented}
  
  \item{clipMode}{\code{character(1)} indicating the clipMode of the returned object.}
  
  \item{verbose}{\code{logical(1)} indicating whether to return verbose output.}


}

\value{A \code{\linkS4class{SffReads}} or \code{\linkS4class{SffReadsQ}} object.}

\author{Matt Settles <msettles@uidaho.edu>}

\examples{
sff <- readSff(system.file("extdata","Small454Test.sff",package="rSFFreader"))
}